#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "analysis.h"

#define PS_PER_TIME_UNIT 1.0967      // one internal time unit in ps (t0 from A1)
#define SITES 5                      // united-atom sites per pentane molecule
#define TRANS_LIMIT (2.0 * PI / 3.0) // |phi| above 120 degrees is trans (barrier top)

// IUPAC dihedral angle, using the minimum-image bond vectors.
double dihedral_angle(Vec3D ri, Vec3D rj, Vec3D rk, Vec3D rl, Vec3D L)
{
    Vec3D b1 = v3_min_image(v3_sub(rj, ri), L);
    Vec3D b2 = v3_min_image(v3_sub(rk, rj), L);
    Vec3D b3 = v3_min_image(v3_sub(rl, rk), L);

    Vec3D n1 = v3_cross(b1, b2);
    Vec3D n2 = v3_cross(b2, b3);

    double y = v3_norm(b2) * v3_dot(b1, n2);
    double x = v3_dot(n1, n2);

    return atan2(y, x);
}

// Mass-weighted centre of mass of every molecule, from the unwrapped positions.
static void molecule_com(struct Parameters *p, struct Vectors *v, Vec3D *com, size_t num_mol)
{
    for (size_t m = 0; m < num_mol; m++)
    {
        Vec3D sum = v3(0.0, 0.0, 0.0);
        double mass_tot = 0.0;

        for (size_t s = 0; s < SITES; s++)
        {
            size_t i = SITES * m + s;
            double mass = p->mass[v->type[i]];
            sum = v3_add(sum, v3_scl(mass, v->r_unwrapped[i]));
            mass_tot += mass;
        }
        com[m] = v3_scl(1.0 / mass_tot, sum);
    }
}

void analysis_init(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    a->num_mol = p->num_part / SITES;
    a->step = 0;
    a->num_blocks = p->num_dt_steps / p->num_dt_block + 1;

    a->hist = (size_t *)calloc(a->num_blocks * PHI_BINS * PHI_BINS, sizeof(size_t));
    a->prev_state = (int *)malloc(2 * a->num_mol * sizeof(int));
    for (size_t k = 0; k < 2 * a->num_mol; k++)
        a->prev_state[k] = -1;
    a->samples_in[0] = a->samples_in[1] = 0;
    a->exits[0] = a->exits[1] = 0;
    a->trace_len = 0;
    a->trace = (double *)malloc(3 * TRACE_MAX * sizeof(double));

    for (size_t i = 0; i < p->num_part; i++)
        v->r_unwrapped[i] = v->r[i];

    a->num_origins = 100;
    a->num_lags = 0;
    a->origins_started = 0;
    a->origin_com = NULL;
    a->origin_step = NULL;
    a->msd_sum = NULL;
    a->msd_cnt = NULL;
    a->msd_block_sum = NULL;
    a->msd_block_cnt = NULL;
    a->com = (Vec3D *)malloc(a->num_mol * sizeof(Vec3D));

    if (p->msd_on)
    {
        // lags 0 .. num_origins * origin interval (in samples) - 1
        a->num_lags = a->num_origins * p->num_dt_msd_origin / p->num_dt_msd;
        a->origin_com = (Vec3D *)malloc(a->num_origins * a->num_mol * sizeof(Vec3D));
        a->origin_step = (size_t *)calloc(a->num_origins, sizeof(size_t));
        a->msd_sum = (double *)calloc(a->num_lags, sizeof(double));
        a->msd_cnt = (size_t *)calloc(a->num_lags, sizeof(size_t));
        a->msd_block_sum = (double *)calloc(a->num_blocks * a->num_lags, sizeof(double));
        a->msd_block_cnt = (size_t *)calloc(a->num_blocks * a->num_lags, sizeof(size_t));
    }
}

// Classify both dihedrals of every molecule: joint histogram, residence counters, trace of molecule 0.
static void sample_dihedrals(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    size_t block = a->step / p->num_dt_block;
    size_t *hist = a->hist + block * PHI_BINS * PHI_BINS;
    Vec3D *r = v->r;
    double phi[2];

    for (size_t m = 0; m < a->num_mol; m++)
    {
        size_t i0 = SITES * m;
        phi[0] = dihedral_angle(r[i0 + 0], r[i0 + 1], r[i0 + 2], r[i0 + 3], p->L);
        phi[1] = dihedral_angle(r[i0 + 1], r[i0 + 2], r[i0 + 3], r[i0 + 4], p->L);

        int bin[2];
        for (int d = 0; d < 2; d++)
        {
            bin[d] = (int)floor((phi[d] + PI) / (2.0 * PI) * PHI_BINS);
            if (bin[d] < 0) bin[d] = 0;
            if (bin[d] >= PHI_BINS) bin[d] = PHI_BINS - 1;

            int state = (fabs(phi[d]) > TRANS_LIMIT) ? 0 : 1;
            int *prev = &a->prev_state[2 * m + d];
            if (*prev >= 0 && *prev != state)
                a->exits[*prev]++;
            *prev = state;
            a->samples_in[state]++;
        }
        hist[bin[0] * PHI_BINS + bin[1]]++;

        if (m == 0 && a->trace_len < TRACE_MAX)
        {
            double *t = a->trace + 3 * a->trace_len++;
            t[0] = a->step * p->dt * PS_PER_TIME_UNIT;
            t[1] = phi[0] * 180.0 / PI;
            t[2] = phi[1] * 180.0 / PI;
        }
    }
}

// Add the squared displacement of the molecular centres of mass from every stored time origin.
static void sample_msd(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    molecule_com(p, v, a->com, a->num_mol);

    // A new time origin every num_dt_msd_origin steps, overwriting the oldest one
    if (a->step % p->num_dt_msd_origin == 0)
    {
        size_t slot = a->origins_started % a->num_origins;
        for (size_t m = 0; m < a->num_mol; m++)
            a->origin_com[slot * a->num_mol + m] = a->com[m];
        a->origin_step[slot] = a->step;
        a->origins_started++;
    }

    size_t active = a->origins_started < a->num_origins ? a->origins_started : a->num_origins;
    for (size_t s = 0; s < active; s++)
    {
        size_t lag = (a->step - a->origin_step[s]) / p->num_dt_msd;
        if (lag >= a->num_lags)
            continue;

        double sum = 0.0;
        for (size_t m = 0; m < a->num_mol; m++)
        {
            Vec3D d = v3_sub(a->com[m], a->origin_com[s * a->num_mol + m]);
            sum += v3_dot(d, d);
        }
        sum /= (double)a->num_mol;

        size_t block = a->origin_step[s] / p->num_dt_block;
        a->msd_sum[lag] += sum;
        a->msd_cnt[lag]++;
        a->msd_block_sum[block * a->num_lags + lag] += sum;
        a->msd_block_cnt[block * a->num_lags + lag]++;
    }
}

void analysis_update(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    // dr holds the displacement of the last step, before the periodic wrap
    for (size_t i = 0; i < p->num_part; i++)
        v->r_unwrapped[i] = v3_add(v->r_unwrapped[i], v->dr[i]);

    a->step++;

    if (a->step % p->num_dt_phi == 0)
        sample_dihedrals(p, v, a);

    if (p->msd_on && a->step % p->num_dt_msd == 0)
        sample_msd(p, v, a);

    if (a->step % p->num_dt_block == 0)
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

void analysis_write(struct Parameters *p, struct Analysis *a)
{
    FILE *fp;
    double ps_per_block = p->num_dt_block * p->dt * PS_PER_TIME_UNIT;
    double ps_per_sample = p->num_dt_phi * p->dt * PS_PER_TIME_UNIT;

    // Blocks that were completed or started
    size_t last = (a->step + p->num_dt_block - 1) / p->num_dt_block;
    if (last > a->num_blocks) last = a->num_blocks;

    if ((fp = open_file(p, "_hist2d.csv")))
    {
        fprintf(fp, "block,t_start_ps,t_end_ps,bin1,bin2,count\n");
        for (size_t b = 0; b < last; b++)
            for (int i = 0; i < PHI_BINS; i++)
                for (int j = 0; j < PHI_BINS; j++)
                {
                    size_t c = a->hist[(b * PHI_BINS + i) * PHI_BINS + j];
                    if (c > 0)
                        fprintf(fp, "%zu,%.4f,%.4f,%d,%d,%zu\n", b, b * ps_per_block, (b + 1) * ps_per_block, i, j, c);
                }
        fclose(fp);
    }

    if ((fp = open_file(p, "_dwell.csv")))
    {
        fprintf(fp, "state,samples,exits,sample_interval_ps\n");
        fprintf(fp, "trans,%zu,%zu,%.6f\n", a->samples_in[0], a->exits[0], ps_per_sample);
        fprintf(fp, "gauche,%zu,%zu,%.6f\n", a->samples_in[1], a->exits[1], ps_per_sample);
        fclose(fp);
    }

    if ((fp = open_file(p, "_phi_trace.csv")))
    {
        fprintf(fp, "time_ps,phi1_deg,phi2_deg\n");
        for (size_t k = 0; k < a->trace_len; k++)
            fprintf(fp, "%.4f,%.3f,%.3f\n", a->trace[3 * k], a->trace[3 * k + 1], a->trace[3 * k + 2]);
        fclose(fp);
    }

    if (p->msd_on && (fp = open_file(p, "_msd.csv")))
    {
        fprintf(fp, "lag_ps,msd_A2,count");
        for (size_t b = 0; b < last; b++)
            fprintf(fp, ",block%zu", b);
        fprintf(fp, "\n");
        for (size_t l = 0; l < a->num_lags; l++)
        {
            if (a->msd_cnt[l] == 0)
                continue;
            fprintf(fp, "%.6f,%.8g,%zu", l * p->num_dt_msd * p->dt * PS_PER_TIME_UNIT,
                    a->msd_sum[l] / (double)a->msd_cnt[l], a->msd_cnt[l]);
            for (size_t b = 0; b < last; b++)
            {
                size_t c = a->msd_block_cnt[b * a->num_lags + l];
                if (c > 0)
                    fprintf(fp, ",%.8g", a->msd_block_sum[b * a->num_lags + l] / (double)c);
                else
                    fprintf(fp, ",");
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

void analysis_free(struct Analysis *a)
{
    free(a->hist);
    free(a->prev_state);
    free(a->trace);
    free(a->origin_com);
    free(a->origin_step);
    free(a->msd_sum);
    free(a->msd_cnt);
    free(a->msd_block_sum);
    free(a->msd_block_cnt);
    free(a->com);
}
