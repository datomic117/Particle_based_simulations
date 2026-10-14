#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "analysis.h"

static FILE *open_file(struct Parameters *p, const char *suffix);
static FILE *open_file_append(struct Parameters *p, const char *suffix);

void analysis_init(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    (void)p;
    (void)v;
    a->step = 0;
    a->rdf_samples = 0;
    a->velocity_frames = 0;
    a->chain_size_samples = 0;
    a->chain_end_to_end_sq_sum = 0.0;
    a->composition_profile_bins = 0;
    a->composition_block_index = 1;
    a->composition_block_samples = 0;
    a->composition_block_start_step = 0;
    a->composition_block_end_step = 0;
    a->composition_count_a = NULL;
    a->composition_count_b = NULL;
    for (size_t t = 0; t < NUM_TYPES; t++)
    {
        a->type_count[t] = 0;
        a->velocity_samples[t] = 0;
        a->velocity_underflow[t] = 0;
        a->velocity_overflow[t] = 0;
        for (size_t b = 0; b < VELOCITY_HIST_BINS; b++)
            a->velocity_hist[t][b] = 0;
        for (size_t component = 0; component < 3; component++)
        {
            a->velocity_component_count[t][component] = 0;
            a->velocity_component_mean[t][component] = 0.0;
            a->velocity_component_m2[t][component] = 0.0;
        }
    }

    for (size_t ti = 0; ti < NUM_TYPES; ti++)
        for (size_t tj = 0; tj < NUM_TYPES; tj++)
            for (size_t b = 0; b < RDF_BINS; b++)
                a->rdf_counts[ti][tj][b] = 0;

    if (p->composition_profile_on)
    {
        if (p->composition_bin_width <= 0.0
            || p->num_dt_composition == 0
            || p->num_dt_composition_block == 0
            || p->L.x <= 0.0
            || p->L.y <= 0.0
            || p->L.z <= 0.0)
        {
            fprintf(stderr, "cannot initialize composition profile: box, bin width, "
                            "and sampling/block intervals must be positive\n");
            exit(EXIT_FAILURE);
        }

        a->composition_profile_bins =
            (size_t)ceil(p->L.x / p->composition_bin_width);
        a->composition_count_a =
            calloc(a->composition_profile_bins, sizeof(*a->composition_count_a));
        a->composition_count_b =
            calloc(a->composition_profile_bins, sizeof(*a->composition_count_b));
        if (a->composition_count_a == NULL || a->composition_count_b == NULL)
        {
            fprintf(stderr, "cannot allocate composition profile bins\n");
            exit(EXIT_FAILURE);
        }

        char filename[1100];
        int filename_length = snprintf(filename, sizeof(filename), "%s%s",
                                       p->filename_analysis,
                                       "_composition_blocks.csv");
        if (filename_length < 0 || (size_t)filename_length >= sizeof(filename))
        {
            fprintf(stderr, "composition profile output filename is too long\n");
            exit(EXIT_FAILURE);
        }

        FILE *fp = fopen(filename, "w");
        if (fp == NULL)
        {
            fprintf(stderr, "cannot open composition profile output '%s'\n", filename);
            exit(EXIT_FAILURE);
        }
        fprintf(fp, "block,step_start,step_end,samples,x,rho_A,rho_B,phi_A\n");
        if (fclose(fp) != 0)
        {
            fprintf(stderr, "cannot finish initializing composition profile output '%s'\n",
                    filename);
            exit(EXIT_FAILURE);
        }
    }
}

static void sample_rdf(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    for (size_t t = 0; t < NUM_TYPES; t++)
        a->type_count[t] = 0;

    for (size_t i = 0; i < p->num_part; i++)
        a->type_count[v->type[i]]++;

    for (int ti = 0; ti < NUM_TYPES; ti++)
    {
        for (int tj = ti; tj < NUM_TYPES; tj++)
        {
            for (size_t i = 0; i < p->num_part; i++)
            {
                if (v->type[i] != ti)
                    continue;

                for (size_t j = (ti == tj) ? i + 1 : 0; j < p->num_part; j++)
                {
                    if (v->type[j] != tj)
                        continue;
                    if (ti == tj && j <= i)
                        continue;

                    Vec3D rij = v3_min_image(v3_sub(v->r[i], v->r[j]), p->L);
                    double r = v3_norm(rij);
                    if (r >= RDF_MAX_R)
                        continue;

                    int bin = (int)(r / RDF_DR);
                    if (bin >= RDF_BINS)
                        bin = RDF_BINS - 1;
                    a->rdf_counts[ti][tj][bin]++;
                }
            }
        }
    }

    a->rdf_samples++;
}

static void sample_velocity(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    const double bin_width = 2.0 * VELOCITY_HIST_RANGE / VELOCITY_HIST_BINS;

    for (size_t i = 0; i < p->num_part; i++)
    {
        int type = v->type[i];
        if (type < 0 || type >= NUM_TYPES)
        {
            fprintf(stderr, "cannot sample velocity: particle %zu has invalid type %d\n", i, type);
            continue;
        }
        if (p->mass[type] <= 0.0 || p->kT <= 0.0)
        {
            fprintf(stderr, "cannot sample velocity: mass and target kT must be positive\n");
            return;
        }

        double thermal_speed = sqrt(p->kT / p->mass[type]);
        double components[3] = {v->v[i].x, v->v[i].y, v->v[i].z};
        for (size_t component = 0; component < 3; component++)
        {
            double normalized_velocity = components[component] / thermal_speed;
            size_t n = ++a->velocity_component_count[type][component];
            double delta = normalized_velocity - a->velocity_component_mean[type][component];
            a->velocity_component_mean[type][component] += delta / (double)n;
            a->velocity_component_m2[type][component] +=
                delta * (normalized_velocity - a->velocity_component_mean[type][component]);

            a->velocity_samples[type]++;
            if (normalized_velocity < -VELOCITY_HIST_RANGE)
            {
                a->velocity_underflow[type]++;
            }
            else if (normalized_velocity >= VELOCITY_HIST_RANGE)
            {
                a->velocity_overflow[type]++;
            }
            else
            {
                size_t bin = (size_t)((normalized_velocity + VELOCITY_HIST_RANGE) / bin_width);
                a->velocity_hist[type][bin]++;
            }
        }
    }

    a->velocity_frames++;
}

static void sample_chain_size(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    if (p->chain_length < 2.0 || floor(p->chain_length) != p->chain_length)
    {
        fprintf(stderr, "cannot sample chain size: chain_length must be an integer >= 2\n");
        return;
    }

    size_t chain_length = (size_t)p->chain_length;
    if (p->num_part % chain_length != 0)
    {
        fprintf(stderr, "cannot sample chain size: num_part must be divisible by chain_length\n");
        return;
    }

    size_t num_chains = p->num_part / chain_length;
    for (size_t chain = 0; chain < num_chains; chain++)
    {
        size_t first = chain * chain_length;
        Vec3D end_to_end = {0.0, 0.0, 0.0};
        for (size_t bond = 0; bond + 1 < chain_length; bond++)
        {
            Vec3D bond_vector = v3_min_image(
                v3_sub(v->r[first + bond + 1], v->r[first + bond]),
                p->L
            );
            end_to_end = v3_add(end_to_end, bond_vector);
        }
        a->chain_end_to_end_sq_sum += v3_dot(end_to_end, end_to_end);
    }

    a->chain_size_samples += num_chains;
}

static void sample_composition_profile(struct Parameters *p,
                                       struct Vectors *v,
                                       struct Analysis *a)
{
    if (a->composition_block_samples == 0)
        a->composition_block_start_step = a->step;
    a->composition_block_end_step = a->step;

    for (size_t i = 0; i < p->num_part; i++)
    {
        int type = v->type[i];
        if (type < 0 || type >= NUM_TYPES)
        {
            fprintf(stderr, "cannot sample composition profile: particle %zu has invalid type %d\n",
                    i, type);
            exit(EXIT_FAILURE);
        }

        double x = v->r[i].x;
        if (x < 0.0 || x >= p->L.x)
            x -= p->L.x * floor(x / p->L.x);
        size_t bin = (size_t)(x / p->composition_bin_width);
        if (bin >= a->composition_profile_bins)
            bin = a->composition_profile_bins - 1;

        if (type == 0)
            a->composition_count_a[bin]++;
        else
            a->composition_count_b[bin]++;
    }

    a->composition_block_samples++;
}

static void write_composition_profile_block(struct Parameters *p, struct Analysis *a)
{
    if (a->composition_block_samples == 0)
        return;

    FILE *fp = open_file_append(p, "_composition_blocks.csv");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    const double cross_section = p->L.y * p->L.z;
    for (size_t bin = 0; bin < a->composition_profile_bins; bin++)
    {
        double x_lo = (double)bin * p->composition_bin_width;
        double x_hi = fmin(x_lo + p->composition_bin_width, p->L.x);
        double x_center = 0.5 * (x_lo + x_hi);
        double slab_volume = (x_hi - x_lo) * cross_section;
        double normalization = (double)a->composition_block_samples * slab_volume;
        double rho_a = (double)a->composition_count_a[bin] / normalization;
        double rho_b = (double)a->composition_count_b[bin] / normalization;
        double density = rho_a + rho_b;
        double phi_a = density > 0.0 ? rho_a / density : NAN;

        fprintf(fp, "%zu,%zu,%zu,%zu,%.8g,%.12g,%.12g,%.12g\n",
                a->composition_block_index,
                a->composition_block_start_step,
                a->composition_block_end_step,
                a->composition_block_samples,
                x_center, rho_a, rho_b, phi_a);
    }

    int write_failed = ferror(fp);
    if (fclose(fp) != 0)
        write_failed = 1;
    if (write_failed)
    {
        fprintf(stderr, "failed writing composition profile block %zu\n",
                a->composition_block_index);
        exit(EXIT_FAILURE);
    }

    a->composition_block_index++;
    a->composition_block_samples = 0;
    a->composition_block_start_step = 0;
    a->composition_block_end_step = 0;
    for (size_t bin = 0; bin < a->composition_profile_bins; bin++)
    {
        a->composition_count_a[bin] = 0;
        a->composition_count_b[bin] = 0;
    }
}

void analysis_update(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    a->step++;

    if (p->analysis_on && p->num_dt_phi > 0 && a->step % p->num_dt_phi == 0)
        sample_rdf(p, v, a);

    if (p->velocity_analysis_on
        && p->num_dt_velocity > 0
        && a->step >= p->velocity_start_step
        && (a->step - p->velocity_start_step) % p->num_dt_velocity == 0)
        sample_velocity(p, v, a);

    if (p->chain_size_analysis_on
        && p->num_dt_chain_size > 0
        && a->step >= p->chain_size_start_step
        && (a->step - p->chain_size_start_step) % p->num_dt_chain_size == 0)
        sample_chain_size(p, v, a);

    if (p->composition_profile_on
        && p->num_dt_composition > 0
        && a->step % p->num_dt_composition == 0)
    {
        sample_composition_profile(p, v, a);
        if (a->step % p->num_dt_composition_block == 0)
            write_composition_profile_block(p, a);
    }

    if (p->num_dt_block > 0 && a->step % p->num_dt_block == 0)
        analysis_write(p, a);
}

static FILE *open_file(struct Parameters *p, const char *suffix)
{
    char name[1100];
    snprintf(name, sizeof(name), "%s%s", p->filename_analysis, suffix);
    FILE *fp = fopen(name, "w");
    if (fp == NULL)
        fprintf(stderr, "cannot open analysis file '%s'\n", name);
    return fp;
}

static FILE *open_file_append(struct Parameters *p, const char *suffix)
{
    char name[1100];
    snprintf(name, sizeof(name), "%s%s", p->filename_analysis, suffix);
    FILE *fp = fopen(name, "a");
    if (fp == NULL)
        fprintf(stderr, "cannot open analysis file '%s' for appending\n", name);
    return fp;
}

static double rdf_shell_volume(double r_lo, double r_hi)
{
    if (r_hi <= r_lo)
        return 0.0;

    r_hi = fmin(r_hi, RDF_MAX_R);

    return (4.0 * PI / 3.0) * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
}

static double expected_pair_count(struct Analysis *a, struct Parameters *p, int ti, int tj, double r_lo, double r_hi)
{
    double volume = p->L.x * p->L.y * p->L.z;
    double shell = rdf_shell_volume(r_lo, r_hi);
    double ni = (double)a->type_count[ti];
    double nj = (double)a->type_count[tj];

    if (ti == tj)
        return 0.5 * ni * (ni - 1.0) * shell / volume;

    return ni * nj * shell / volume;
}

static const char *particle_type_name(int type)
{
    static const char *const names[] = {"A", "B"};
    return type < 2 ? names[type] : "unknown";
}

static void write_velocity_histogram(struct Parameters *p, struct Analysis *a)
{
    FILE *fp = open_file(p, "_velocity_hist.csv");
    if (!fp)
        return;

    fprintf(fp, "type,velocity,density,maxwell_boltzmann,count,samples,underflow,overflow\n");
    const double bin_width = 2.0 * VELOCITY_HIST_RANGE / VELOCITY_HIST_BINS;
    for (int type = 0; type < NUM_TYPES; type++)
    {
        if (a->velocity_samples[type] == 0 || p->mass[type] <= 0.0 || p->kT <= 0.0)
            continue;

        double thermal_speed = sqrt(p->kT / p->mass[type]);
        double velocity_bin_width = bin_width * thermal_speed;
        for (size_t bin = 0; bin < VELOCITY_HIST_BINS; bin++)
        {
            double normalized_velocity = -VELOCITY_HIST_RANGE + ((double)bin + 0.5) * bin_width;
            double velocity = normalized_velocity * thermal_speed;
            double density = (double)a->velocity_hist[type][bin]
                / ((double)a->velocity_samples[type] * velocity_bin_width);
            double maxwell_boltzmann = sqrt(p->mass[type] / (2.0 * PI * p->kT))
                * exp(-p->mass[type] * velocity * velocity / (2.0 * p->kT));
            fprintf(fp, "%s,%.12g,%.12g,%.12g,%zu,%zu,%zu,%zu\n",
                    particle_type_name(type), velocity, density, maxwell_boltzmann,
                    a->velocity_hist[type][bin], a->velocity_samples[type],
                    a->velocity_underflow[type], a->velocity_overflow[type]);
        }
    }
    fclose(fp);

    fp = open_file(p, "_velocity_stats.csv");
    if (!fp)
        return;

    fprintf(fp, "type,component,samples,mean_velocity,variance,kBT_measured,kBT_target\n");
    static const char *const component_names[] = {"x", "y", "z"};
    for (int type = 0; type < NUM_TYPES; type++)
    {
        for (size_t component = 0; component < 3; component++)
        {
            size_t n = a->velocity_component_count[type][component];
            if (n < 2 || p->mass[type] <= 0.0 || p->kT <= 0.0)
                continue;

            double thermal_speed = sqrt(p->kT / p->mass[type]);
            double variance_normalized =
                a->velocity_component_m2[type][component] / (double)(n - 1);
            double variance = variance_normalized * thermal_speed * thermal_speed;
            double measured_kT = p->mass[type] * variance;
            double mean_velocity =
                a->velocity_component_mean[type][component] * thermal_speed;
            fprintf(fp, "%s,%s,%zu,%.12g,%.12g,%.12g,%.12g\n",
                    particle_type_name(type), component_names[component], n,
                    mean_velocity, variance, measured_kT, p->kT);
        }
    }
    fclose(fp);
}

void analysis_write(struct Parameters *p, struct Analysis *a)
{
    if (p->composition_profile_on)
        write_composition_profile_block(p, a);

    if (a->chain_size_samples > 0)
    {
        FILE *fp = open_file(p, "_chain_size.csv");
        if (fp)
        {
            size_t chain_length = (size_t)p->chain_length;
            double measured = a->chain_end_to_end_sq_sum / (double)a->chain_size_samples;
            double ideal = 3.0 * (double)(chain_length - 1) * p->kT / p->k_b;
            fprintf(fp, "chains,chain_length,samples,mean_square_end_to_end,ideal_chain_mean_square\n");
            fprintf(fp, "%zu,%zu,%zu,%.12g,%.12g\n",
                    p->num_part / chain_length, chain_length, a->chain_size_samples,
                    measured, ideal);
            fclose(fp);
        }
    }

    if (a->rdf_samples > 0)
    {
        FILE *fp = open_file(p, "_rdf.csv");
        if (fp)
        {
            fprintf(fp, "pair,r_min,r_max,g_r,count,expected\n");

            for (int ti = 0; ti < NUM_TYPES; ti++)
            {
                for (int tj = ti; tj < NUM_TYPES; tj++)
                {
                    const char *pair_name = "pair";
                    if (NUM_TYPES == 2)
                    {
                        if (ti == 0 && tj == 0) pair_name = "AA";
                        else if (ti == 0 && tj == 1) pair_name = "AB";
                        else if (ti == 1 && tj == 1) pair_name = "BB";
                    }

                    for (int b = 0; b < RDF_BINS; b++)
                    {
                        double r_lo = (double)b * RDF_DR;
                        double r_hi = r_lo + RDF_DR;

                        double expected = expected_pair_count(a, p, ti, tj, r_lo, r_hi);
                        double count = (double)a->rdf_counts[ti][tj][b];
                        double g = expected > 0.0 ? count / ((double)a->rdf_samples * expected) : 0.0;

                        fprintf(fp, "%s,%.3f,%.3f,%.8f,%zu,%.8f\n",
                                pair_name, r_lo, r_hi, g, a->rdf_counts[ti][tj][b], expected);
                    }
                }
            }
            fclose(fp);
        }
    }

    if (a->velocity_frames > 0)
        write_velocity_histogram(p, a);
}

void analysis_free(struct Analysis *a)
{
    free(a->composition_count_a);
    free(a->composition_count_b);
    a->composition_count_a = NULL;
    a->composition_count_b = NULL;
}
