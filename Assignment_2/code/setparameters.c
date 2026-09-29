#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// B3: continue equilibration for 50 ps from the latest saved state.
void set_parameters(struct Parameters *p_parameters)
{
    // Target temperature: 293 K.
    p_parameters->kT = 293.0;

    // Thermostat coupling time: 0.1 ps in internal units.
    p_parameters->tau_T = 0.0911836;

    // Masses in atomic mass units.
    p_parameters->mass[TYPE_CH3] = 15.035;
    p_parameters->mass[TYPE_CH2] = 14.027;

    // Lennard-Jones energy parameters: epsilon / k_B in Kelvin.
    p_parameters->epsilon[TYPE_CH3] = 98.0;
    p_parameters->epsilon[TYPE_CH2] = 46.0;

    // Lennard-Jones length parameters in angstroms.
    p_parameters->sigma[TYPE_CH3] = 3.75;
    p_parameters->sigma[TYPE_CH2] = 3.95;

    // Particle count must match the saved configuration.
    p_parameters->num_part = 2000;

    // Continue the simulation with the thermostat enabled.
    p_parameters->force_test = 0;
    p_parameters->load_restart = 1;
    p_parameters->is_NVT = 1;

    // 50,000 steps at 1 fs = 50 ps.
    p_parameters->num_dt_steps = 50000;
    p_parameters->num_dt_output = 100;

    // Exclusions used once molecular connectivity is implemented.
    p_parameters->factor_12_nb = 0.0;
    p_parameters->factor_13_nb = 0.0;
    p_parameters->factor_14_nb = 0.0;

    // Timestep: 1 fs in internal units.
    p_parameters->dt = 0.000911836;

    // Box size at a mass density of 626 kg/m^3.
    const double molecule_mass =
        2.0 * p_parameters->mass[TYPE_CH3]
        + 3.0 * p_parameters->mass[TYPE_CH2];

    const double total_mass_kg =
        (p_parameters->num_part / 5)
        * molecule_mass * 1.66053906660e-27;

    const double box_length_angstrom =
        cbrt(total_mass_kg / 626.0) / 1.0e-10;

    p_parameters->L = (struct Vec3D){
        box_length_angstrom,
        box_length_angstrom,
        box_length_angstrom
    };

    // Interaction cutoff and neighbour-list buffer in angstroms.
    p_parameters->r_cut = 14.0;
    p_parameters->r_shell = 2.0;

    // Trajectory output.
    p_parameters->num_dt_pdb = 500;
    p_parameters->rescale_output = 1.0;

    strcpy(p_parameters->filename_pdb, "../data/b3_equil_long");
    strcpy(p_parameters->filename_xyz, "../data/b3_equil_long");

    // Continue from the completed extra equilibration run.
    strcpy(
        p_parameters->restart_in_filename,
        "../data/b3_equil_extra_restart.dat"
    );

    // Save to a separate file to preserve the starting configuration.
    p_parameters->num_dt_restart = 1000;

    strcpy(
        p_parameters->restart_out_filename,
        "../data/b3_equil_long_restart.dat"
    );

    // Minimum-image requirement.
    if (p_parameters->r_cut > p_parameters->L.x / 2.0 ||
        p_parameters->r_cut > p_parameters->L.y / 2.0 ||
        p_parameters->r_cut > p_parameters->L.z / 2.0)
    {
        fprintf(stderr, "Error: cutoff exceeds half the box length.\n");
        exit(EXIT_FAILURE);
    }
}