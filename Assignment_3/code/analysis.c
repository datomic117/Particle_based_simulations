#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "analysis.h"

void analysis_init(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    (void)p;
    (void)v;
    a->step = 0;
    a->rdf_samples = 0;
    a->rdf_block_samples = 0;
    a->rdf_blocks_written = 0;
    a->profile_samples = 0;
    a->profile_blocks_written = 0;
    a->fp_pressure = NULL;
    for (size_t t = 0; t < NUM_TYPES; t++)
    {
        for (size_t b = 0; b < PROFILE_BINS_MAX; b++)
            a->profile_counts[t][b] = 0;
        for (size_t u = 0; u < NUM_TYPES; u++)
            for (size_t b = 0; b < RDF_BINS; b++)
                a->rdf_block_counts[t][u][b] = 0;
    }
    a->velocity_frames = 0;
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
                    a->rdf_block_counts[ti][tj][bin]++;
                }
            }
        }
    }

    a->rdf_samples++;
    a->rdf_block_samples++;
}

// Number of x-bins of the composition profile and the bin width actually used
static size_t profile_bins(struct Parameters *p, double *dx)
{
    size_t n = (size_t)(p->L.x / PROFILE_DX + 0.5);
    if (n > PROFILE_BINS_MAX)
        n = PROFILE_BINS_MAX;
    *dx = p->L.x / (double)n;
    return n;
}

// Add the current positions to the composition profile of this block
static void sample_profile(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    double dx;
    size_t n = profile_bins(p, &dx);

    for (size_t i = 0; i < p->num_part; i++)
    {
        size_t bin = (size_t)(v->r[i].x / dx);
        if (bin >= n)
            bin = n - 1;
        a->profile_counts[v->type[i]][bin]++;
    }
    a->profile_samples++;
}

// Append the diagonal pressure tensor of this step: P_aa = (sum m v_a^2 + sum r_a F_a)/V
static void sample_pressure(struct Parameters *p, struct Vectors *v, struct Analysis *a, double time)
{
    if (a->fp_pressure == NULL)
    {
        char name[1100];
        snprintf(name, sizeof(name), "%s_pressure.csv", p->filename_analysis);
        a->fp_pressure = fopen(name, "w");
        if (a->fp_pressure == NULL)
        {
            fprintf(stderr, "cannot open analysis file '%s'\n", name);
            return;
        }
        fprintf(a->fp_pressure, "step,time,Pxx,Pyy,Pzz\n");
    }
    fprintf(a->fp_pressure, "%zu,%.8g,%.10g,%.10g,%.10g\n", a->step, time,
            v->press_kin_xyz.x + v->press_vir_xyz.x,
            v->press_kin_xyz.y + v->press_vir_xyz.y,
            v->press_kin_xyz.z + v->press_vir_xyz.z);
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

static void write_profile_block(struct Parameters *p, struct Analysis *a);
static void write_rdf_block(struct Parameters *p, struct Analysis *a);

void analysis_update(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    a->step++;

    if (p->analysis_on && p->num_dt_phi > 0 && a->step >= p->analysis_start_step
        && a->step % p->num_dt_phi == 0)
        sample_rdf(p, v, a);

    if (p->profile_on && p->num_dt_profile > 0 && a->step % p->num_dt_profile == 0)
        sample_profile(p, v, a);

    if (p->profile_on && p->num_dt_pressure > 0 && a->step % p->num_dt_pressure == 0)
        sample_pressure(p, v, a, (double)a->step * p->dt);

    if (p->velocity_analysis_on
        && p->num_dt_velocity > 0
        && a->step >= p->velocity_start_step
        && (a->step - p->velocity_start_step) % p->num_dt_velocity == 0)
        sample_velocity(p, v, a);

    if (p->num_dt_block > 0 && a->step % p->num_dt_block == 0)
    {
        analysis_write(p, a);
        write_profile_block(p, a);
        write_rdf_block(p, a);
    }
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

// Append the block-averaged number densities rho_A(x), rho_B(x) of the block that just ended
static void write_profile_block(struct Parameters *p, struct Analysis *a)
{
    if (!p->profile_on || a->profile_samples == 0)
        return;

    char name[1100];
    snprintf(name, sizeof(name), "%s_profile.csv", p->filename_analysis);
    FILE *fp = fopen(name, a->profile_blocks_written == 0 ? "w" : "a");
    if (fp == NULL)
    {
        fprintf(stderr, "cannot open analysis file '%s'\n", name);
        return;
    }
    if (a->profile_blocks_written == 0)
        fprintf(fp, "block,step_end,x,rho_A,rho_B\n");

    double dx;
    size_t n = profile_bins(p, &dx);
    double slab = dx * p->L.y * p->L.z;
    for (size_t b = 0; b < n; b++)
    {
        fprintf(fp, "%zu,%zu,%.6f,%.6f,%.6f\n", a->profile_blocks_written, a->step, (b + 0.5) * dx,
                (double)a->profile_counts[0][b] / ((double)a->profile_samples * slab),
                (double)a->profile_counts[1][b] / ((double)a->profile_samples * slab));
        a->profile_counts[0][b] = 0;
        a->profile_counts[1][b] = 0;
    }
    fclose(fp);
    a->profile_blocks_written++;
    a->profile_samples = 0;
}

// Append the g(r) of the block that just ended (same normalisation as the cumulative g(r))
static void write_rdf_block(struct Parameters *p, struct Analysis *a)
{
    if (!p->analysis_on || a->rdf_block_samples == 0)
        return;

    char name[1100];
    snprintf(name, sizeof(name), "%s_rdf_blocks.csv", p->filename_analysis);
    FILE *fp = fopen(name, a->rdf_blocks_written == 0 ? "w" : "a");
    if (fp == NULL)
    {
        fprintf(stderr, "cannot open analysis file '%s'\n", name);
        return;
    }
    if (a->rdf_blocks_written == 0)
        fprintf(fp, "block,pair,r_min,r_max,g_r\n");

    for (int ti = 0; ti < NUM_TYPES; ti++)
    {
        for (int tj = ti; tj < NUM_TYPES; tj++)
        {
            const char *pair_name = ti == 0 ? (tj == 0 ? "AA" : "AB") : "BB";
            for (int b = 0; b < RDF_BINS; b++)
            {
                double r_lo = (double)b * RDF_DR;
                double expected = expected_pair_count(a, p, ti, tj, r_lo, r_lo + RDF_DR);
                double g = expected > 0.0
                    ? (double)a->rdf_block_counts[ti][tj][b] / ((double)a->rdf_block_samples * expected)
                    : 0.0;
                fprintf(fp, "%zu,%s,%.3f,%.3f,%.8f\n", a->rdf_blocks_written, pair_name, r_lo, r_lo + RDF_DR, g);
                a->rdf_block_counts[ti][tj][b] = 0;
            }
        }
    }
    fclose(fp);
    a->rdf_blocks_written++;
    a->rdf_block_samples = 0;
}

void analysis_write(struct Parameters *p, struct Analysis *a)
{
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
    if (a->fp_pressure != NULL)
        fclose(a->fp_pressure);
    a->fp_pressure = NULL;
}
