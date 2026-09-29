#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation. This is the only place run settings
// live: there are no input files, so changing a run means editing the values
// below and recompiling.

void set_parameters(struct Parameters *p_parameters)
{
    /// \todo Replace these demonstration values by the parameters of the n-pentane
    /// force field: type-dependent masses and LJ parameters, the bond, angle and
    /// dihedral parameters, and the thermostat relaxation time

    // The values below are demonstration values in reduced Lennard-Jones units

    p_parameters->mass[0] = 15.035;   // CH3
    p_parameters->mass[1] = 14.027;   // CH2

    p_parameters->epsilon[0] = 98.0;  // CH3, in kB K
    p_parameters->epsilon[1] = 46.0;  // CH2, in kB K

    p_parameters->sigma[0] = 3.75;    // CH3, Angstrom
    p_parameters->sigma[1] = 3.95;    // CH2, Angstrom

    // Target thermal energy. Since energy is expressed in kB K,
    // the numerical value equals the temperature in Kelvin.
    p_parameters->kT = 293.0;

    // Berendsen thermostat relaxation time
    p_parameters->tau_T = 0.009118;


    // -------------------------------------------------------------------------
    // Bonded force-field parameters for n-pentane
    // -------------------------------------------------------------------------

    // Harmonic bond potential:
    //
    // U_bond(r) = 1/2 * k_b * (r - r_0)^2
    //
    // Lengths are expressed in Angstrom and energies in kB K, so the numerical
    // values supplied in the assignment can be used directly.
    p_parameters->r_0 = 1.54;          // equilibrium bond length, Angstrom
    p_parameters->k_b = 3.19e5;        // bond force constant, K / Angstrom^2


    // Harmonic angle potential:
    //
    // U_angle(theta) = 1/2 * k_theta * (theta - theta_0)^2
    //
    // Angles must be expressed in radians in the force calculation.
    p_parameters->theta_0 =
        114.0 * M_PI / 180.0;          // equilibrium angle, radians

    p_parameters->k_theta =
        6.25e4;                        // angle force constant, K / rad^2


    // Ryckaert-Bellemans torsion potential:
    //
    // U_tors(phi) =
    //     c_0
    //   + c_1 * cos(phi)
    //   + c_2 * cos(phi)^2
    //   + c_3 * cos(phi)^3
    //
    // The coefficients are given in kB K and can therefore also be used
    // directly in the internal energy units of the program.
    p_parameters->c_0 = 1010.0;
    p_parameters->c_1 = -2018.9;
    p_parameters->c_2 = 136.4;
    p_parameters->c_3 = 3165.3;


    // The parameters below control core functionalities of the code, but many values will need to be changed
    //
    // Run modes:
    //  - production:      force_test = 0, num_dt_steps > 0, is_NVT = 1
    //  - force check:     force_test > 0 (finite-difference test of the forces on every
    //                     force_test-th particle and of the virial, after which the
    //                     program exits; num_dt_steps is ignored)
    //  - energy check:    is_NVT = 0 (NVE: thermostat off, so Etot must be conserved)


    // -------------------------------------------------------------------------
    // B5: long-time NVE energy-conservation test for one pentane molecule
    // -------------------------------------------------------------------------

    // One n-pentane molecule contains five united-atom sites:
    //
    // CH3 - CH2 - CH2 - CH2 - CH3
    p_parameters->num_part = 5;


    // We are no longer doing the finite-difference force test.
    // Run the actual molecular dynamics trajectory.
    p_parameters->force_test = 0;


    // NVE ensemble:
    //
    // no thermostat is applied, so the total energy should remain conserved
    // apart from the bounded numerical error of velocity-Verlet integration.
    p_parameters->is_NVT = 0;


    // Scaling of the non-bonded interaction between nearby atoms of one chain:
    //
    // 1-2, 1-3 and 1-4 interactions are excluded.
    // The 1-5 interaction remains and is therefore included normally.
    p_parameters->factor_12_nb = 0.0;
    p_parameters->factor_13_nb = 0.0;
    p_parameters->factor_14_nb = 0.0;


    /// \todo Set the time step, box size and cut-off distance to values appropriate
    /// for n-pentane at a mass density of 626 kg/m3

    // 1 fs expressed in the internal simulation time unit.
    p_parameters->dt = 0.0009118;


    // Long single-molecule NVE trajectory.
    //
    // 100000 steps at 1 fs = approximately 100 ps.
    p_parameters->num_dt_steps = 100000;


    // Save thermodynamic quantities every 20 steps.
    //
    // This gives 5000 samples over the NVE trajectory, which is more than
    // enough to inspect the total-energy conservation.
    p_parameters->num_dt_output = 20;


    // Use the same large box as for the single-molecule thermalisation.
    //
    // The goal here is an isolated molecule rather than a liquid at the
    // target density.
    p_parameters->L =
        (struct Vec3D){40.0, 40.0, 40.0};


    // Non-bonded cutoff.
    p_parameters->r_cut = 14.0;

    // Neighbor-list shell thickness.
    p_parameters->r_shell = 0.4;


    // A PDB frame is not needed very frequently for the energy test.
    p_parameters->num_dt_pdb = 1000;

    strcpy(
        p_parameters->filename_pdb,
        "b5_nve"
    );                                           // filename (without extension) for pdb file

    p_parameters->rescale_output = 1;            // factor used to rescale output lengthscale
                                                 // (Most visualisation programs identify bonds
                                                 // based on distances of order 1)


    // Start the NVE run from the previously thermalised molecule.
    p_parameters->load_restart = 1;

    strcpy(
        p_parameters->restart_in_filename,
        "restart_b5_thermalised.dat"
    );


    // Save the final NVE configuration separately.
    //
    // Do NOT overwrite restart_b5_thermalised.dat, because that is our common
    // starting point for the B5 verification tests.
    p_parameters->num_dt_restart = 100000;

    strcpy(
        p_parameters->restart_out_filename,
        "restart_b5_nve.dat"
    );


    // The minimum image convention requires r_cut <= L/2: beyond that a particle
    // would interact with two images of the same neighbor.

    if (p_parameters->r_cut > p_parameters->L.x / 2.0)
        fprintf(stderr, "Warning! r_cut > Lx/2");

    if (p_parameters->r_cut > p_parameters->L.y / 2.0)
        fprintf(stderr, "Warning! r_cut > Ly/2");

    if (p_parameters->r_cut > p_parameters->L.z / 2.0)
        fprintf(stderr, "Warning! r_cut > Lz/2");
}