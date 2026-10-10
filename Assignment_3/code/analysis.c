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
    for (size_t t = 0; t < NUM_TYPES; t++)
        a->type_count[t] = 0;

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
                }
            }
        }
    }

    a->rdf_samples++;
}

void analysis_update(struct Parameters *p, struct Vectors *v, struct Analysis *a)
{
    a->step++;

    if (a->step % p->num_dt_phi == 0)
        sample_rdf(p, v, a);

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

    if (a->rdf_samples > 0 && (fp = open_file(p, "_rdf.csv")))
    {
        fprintf(fp, "pair,r_min,r_max,g_r,count,expected\n");
        for (int ti = 0; ti < NUM_TYPES; ti++)
        {
            for (int tj = ti; tj < NUM_TYPES; tj++)
            {
                const char *pair_name = (ti == 0 && tj == 0) ? "AA" : (ti == 0 && tj == 1) ? "AB" : "BB";
                for (int b = 0; b < RDF_BINS; b++)
                {
                    double r_lo = b * RDF_DR;
                    double r_hi = r_lo + RDF_DR;
                    double shell = 4.0 * PI / 3.0 * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
                    double volume = p->L.x * p->L.y * p->L.z;
                    double expected = 0.0;

                    if (ti == tj)
                    {
                        double n_i = (double)a->type_count[ti];
                        expected = 0.5 * n_i * (n_i - 1.0) * shell / volume;
                    }
                    else
                    {
                        expected = (double)a->type_count[ti] * (double)a->type_count[tj] * shell / volume;
                    }

                    double g = expected > 0.0 ? (double)a->rdf_counts[ti][tj][b] / ((double)a->rdf_samples * expected) : 0.0;
                    fprintf(fp, "%s,%.3f,%.3f,%.8f,%zu,%.8f\n",
                            pair_name, r_lo, r_hi, g, a->rdf_counts[ti][tj][b], expected);
                }
            }
        }
        fclose(fp);
    }
}

void analysis_free(struct Analysis *a)
{
    (void)a;
}
