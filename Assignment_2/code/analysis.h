#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <stddef.h>
#include "structs.h"

/// Number of bins per dihedral axis in the joint (phi1, phi2) histogram (5 degrees each).
#define PHI_BINS 72

/// Upper limit on the number of samples stored in the single-molecule phi trace.
#define TRACE_MAX 5000

/**
 * @brief State of the on-the-fly analysis: C1 (dihedral statistics) and C2
 * (mean-square displacement of the molecular centres of mass).
 */
struct Analysis
{
    size_t num_mol;           //!< number of molecules (num_part / 5)
    size_t step;              //!< steps since the analysis was started
    size_t num_blocks;        //!< number of blocks allocated

    // C1: dihedral statistics, per block of num_dt_block steps
    size_t *hist;             //!< [block][PHI_BINS][PHI_BINS] counts of (phi1, phi2)
    int *prev_state;          //!< [2 * num_mol] last state (0 trans, 1 gauche, -1 unknown) of each dihedral
    size_t samples_in[2];     //!< phi samples spent in trans / gauche
    size_t exits[2];          //!< transitions out of trans / gauche
    size_t trace_len;         //!< number of stored trace samples
    double *trace;            //!< [TRACE_MAX][3] time (ps), phi1, phi2 of molecule 0

    // C2: mean-square displacement of the molecular centres of mass
    size_t num_lags;          //!< number of lags stored
    size_t num_origins;       //!< ring size of stored time origins
    size_t origins_started;   //!< number of origins created so far
    Vec3D *origin_com;        //!< [num_origins][num_mol] centres of mass at the origins
    size_t *origin_step;      //!< step of each stored origin
    double *msd_sum;          //!< [num_lags] sum of the molecule-averaged squared displacements
    size_t *msd_cnt;          //!< [num_lags] number of contributions
    double *msd_block_sum;    //!< [block][num_lags], blocks are labelled by the origin's block
    size_t *msd_block_cnt;    //!< [block][num_lags]
    Vec3D *com;               //!< [num_mol] work array: current centres of mass
};

/**
 * @brief Allocate the analysis and start the unwrapped trajectory from the current positions.
 *
 * @param p_parameters used members: num_part, num_dt_steps, num_dt_block, num_dt_msd, num_dt_msd_origin, msd_on
 * @param p_vectors used members: r (copied to r_unwrapped)
 * @param p_analysis analysis state to initialise
 */
void analysis_init(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Analysis *p_analysis);

/**
 * @brief Call once per step after update_positions: accumulates the unwrapped
 * displacement and, at the sampling intervals, updates the dihedral
 * histograms, the residence-time counters and the MSD.
 */
void analysis_update(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Analysis *p_analysis);

/// @brief Write all analysis files (prefix given by Parameters.filename_analysis).
void analysis_write(struct Parameters *p_parameters, struct Analysis *p_analysis);

/// @brief Free the memory held by the analysis.
void analysis_free(struct Analysis *p_analysis);

/**
 * @brief IUPAC dihedral angle of the atoms i-j-k-l in radians in (-pi, pi]; trans = +-pi.
 * Its cosine equals the cos(phi) used by calculate_forces_dihedral.
 */
double dihedral_angle(struct Vec3D ri, struct Vec3D rj, struct Vec3D rk, struct Vec3D rl, struct Vec3D L);

#endif /* ANALYSIS_H_ */
