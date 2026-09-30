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

    // -------------------------------------------------------------------------
    // B9: second Berendsen thermostat test
    // -------------------------------------------------------------------------
    //
    // Weaker thermostat coupling:
    //
    // tau_T = 0.9118 internal time units
    //       ~= 1.0 ps
    //
    // This run will be compared with the strong-coupling case tau = 0.01 ps
    // and with the NVE reference from B7.
    p_parameters->tau_T = 0.9118;


    // -------------------------------------------------------------------------
    // Bonded force-field parameters for n-pentane
    // -------------------------------------------------------------------------

    // Harmonic bond potential:
    //
    // U_bond(r) = 1/2 * k_b * (r - r_0)^2
    p_parameters->r_0 = 1.54;          // Angstrom
    p_parameters->k_b = 3.19e5;        // K / Angstrom^2


    // Harmonic angle potential:
    //
    // U_angle(theta) = 1/2 * k_theta * (theta - theta_0)^2
    p_parameters->theta_0 =
        114.0 * M_PI / 180.0;          // radians

    p_parameters->k_theta =
        6.25e4;                        // K / rad^2


    // Ryckaert-Bellemans torsion potential
    p_parameters->c_0 = 1010.0;
    p_parameters->c_1 = -2018.9;
    p_parameters->c_2 = 136.4;
    p_parameters->c_3 = 3165.3;


    // -------------------------------------------------------------------------
    // B9: NVT test of the full packed pentane system
    // -------------------------------------------------------------------------

    // 400 pentane molecules = 2000 united-atom sites.
    p_parameters->num_part = 2000;


    // Not a force test.
    p_parameters->force_test = 0;


    // NVT ensemble: thermostat enabled.
    p_parameters->is_NVT = 1;


    // Exclude intramolecular 1-2, 1-3 and 1-4 LJ interactions.
    p_parameters->factor_12_nb = 0.0;
    p_parameters->factor_13_nb = 0.0;
    p_parameters->factor_14_nb = 0.0;


    /// \todo Set the time step, box size and cut-off distance to values appropriate
    /// for n-pentane at a mass density of 626 kg/m3

    // 1 fs in internal simulation units.
    p_parameters->dt = 0.0009118;


    // Cubic box corresponding to rho = 626 kg/m^3 for 400 molecules.
    p_parameters->L =
        (struct Vec3D){
            42.4612104,
            42.4612104,
            42.4612104
        };


    // 20000 steps * 1 fs ~= 20 ps.
    p_parameters->num_dt_steps = 20000;


    // Save thermodynamic data every 10 steps.
    p_parameters->num_dt_output = 10;


    // Non-bonded cutoff.
    p_parameters->r_cut = 14.0;

    // Neighbor-list shell.
    p_parameters->r_shell = 0.4;


    // Save PDB frame every 1000 steps.
    p_parameters->num_dt_pdb = 1000;

    strcpy(
        p_parameters->filename_pdb,
        "b9_nvt_tau1"
    );

    p_parameters->rescale_output = 1;


    // Reproduce the same B6 packed starting configuration.
    p_parameters->load_restart = 0;

    strcpy(
        p_parameters->restart_in_filename,
        "restart_b5_thermalised.dat"
    );


    // Save the final state of this second B9 run.
    p_parameters->num_dt_restart = 20000;

    strcpy(
        p_parameters->restart_out_filename,
        "restart_b9_tau1.dat"
    );


    // Minimum-image safety check.
    if (p_parameters->r_cut > p_parameters->L.x / 2.0)
        fprintf(stderr, "Warning! r_cut > Lx/2");

    if (p_parameters->r_cut > p_parameters->L.y / 2.0)
        fprintf(stderr, "Warning! r_cut > Ly/2");

    if (p_parameters->r_cut > p_parameters->L.z / 2.0)
        fprintf(stderr, "Warning! r_cut > Lz/2");
}