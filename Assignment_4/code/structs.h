#ifndef TYPES_MD_H_
#define TYPES_MD_H_
#include <stdbool.h>
#include "constants.h"

/* This header file contains definitions of struct types used in the DEM code */

/**
 * @brief Struct to store i, j, k indices of a 3D grid
 *
 */
struct Index3D
{
    size_t i, j, k; //!< 3 indices: i,j, k
};

/**
 * @brief A 3D vector. All positions, velocities, forces, angular velocities
 * and torques are stored as arrays of Vec3D.
 */
struct Vec3D
{
    double x, y, z; //!< The three Cartesian components of the vector
};

/**
 * @brief A 3D vector stored together with its square length. Used for
 * connecting vectors and accumulated displacements, where the square length is
 * needed anyway for comparisons against a cut-off distance.
 */
struct DeltaR
{
    double x, y, z; //!< x, y and z components
    double sq;      //!< square length x^2 + y^2 + z^2
};

/**
 * @brief All run settings in one struct, set by @ref set_parameters. There are
 * no input files: changing a run means editing setparameters.c and recompiling.
 */
struct Parameters
{
    size_t num_part;       //!< Number of particles
    size_t num_dt_steps;   //!< Number of time steps
    double dt;             //!< integration time step
    struct Vec3D L;        //!< Box size in the three directions
    double Tg;             //!< Granular temperature. Can be used to initialize velocities. 1.5Tg is the average kinetic energy per particle.
    double density;        //!< Density of particles
    double mass_ref;       //!< Reference mass used for collision parameters
    double R_max;          //!< Maximum sphere radius
    double R_min;          //!< Minumum sphere radius 
    double freq;           //!< Frequency of the vertical oscillation of the box, zero for a box at rest. See the vibrated-bed case in cases/DEM_BNE.
    double A;              //!< Amplitude of the vertical oscillation of the box, zero for a box at rest.
    int run;               //!< Which set-up to simulate: RUN_COLLISION, RUN_SLOPE or RUN_SILO (constants.h)
    double v_impact;       //!< Part A: the speed with which a particle hits the wall, or each of two particles the other
    double t_tilt;         //!< Part B: time at which gravity starts to rotate; before it the bed settles under vertical gravity
    double tilt_rate;      //!< Part B: rate (rad/s) at which the slope angle grows from zero after t_tilt
    double slope_angle;    //!< Part B: the current angle (rad) over which gravity is rotated about the y-axis
    double R_cyl;          //!< Parts C, D: radius of the cylindrical silo wall
    double z_silo;         //!< Parts C, D: height of the silo bottom plate; below it is the catch region for the discharged material
    double R_orifice;      //!< Parts C, D: radius of the circular orifice in the silo bottom plate
    int orifice_open;      //!< Parts C, D: 0 while the orifice is closed (filling, settling), 1 once it is open
    double t_open;         //!< Part D: time at which the orifice opens; before it the fill settles on the closed plate
    struct Vec3D g;        //!< gravitational acceleration vector
    double k_n_pp;         //!< normal elastic spring constant for particle-particle interactions
    double eta_n_pp;       //!< normal dashpot friction constant for particle-particle interactions
    double k_t_pp;         //!< tangential elastic spring constant for particle-particle interactions
    double eta_t_pp;       //!< tangential dashpot damping coeff. for particle-particle interactions
    double fric_pp;        //!< friction coefficient for particle-particle interactions
    unsigned int num_walls;//!< number of solid walls in the system (at most NUM_WALLS_MAX)
    bool (*wall_function[NUM_WALLS_MAX])(struct Parameters *, double, struct Vec3D *, struct DeltaR *, struct Vec3D *); //!< the walls, as function pointers; see walls.c for the signature a wall implements and how to register one
    double k_n_pw[NUM_WALLS_MAX];    //!< normal elastic spring constant for particle-wall interactions
    double eta_n_pw[NUM_WALLS_MAX];  //!< normal dashpot damping coeff. for particle-wall interactions
    double k_t_pw[NUM_WALLS_MAX];    //!< tangential elastic spring constant for particle-wall interactions
    double eta_t_pw[NUM_WALLS_MAX];  //!< tangential damping coeff. for particle-wall interactions
    double fric_pw[NUM_WALLS_MAX];   //!< friction coeff. for particle-wall interactions
    double e_n_pp_set;               //!< Prescribed coefficient of normal restitution for particle-particle contacts, the value the contact parameters were tuned to
    double e_n_pw_set[NUM_WALLS_MAX];//!< Prescribed coefficient of normal restitution for particle-wall contacts
    double r_cut;                    //!< Cut-off distance of the pair interaction: the largest centre-to-centre distance at which two particles can touch, 2*R_max
    double r_shell;                  //!< Shell thickness for neighbor list
    size_t num_dt_printf;            //!< Number of time steps between prints to screen
    size_t num_dt_traj;              //!< Number of time steps between trajectory saves
    char filename_xyz[1024];         //!< filename (without extension) for the xyz trajectory file
    char load_restart;               //!< if equal 1 restart file is loaded
    size_t num_dt_restart;           //!< Number of time steps between saves of restart file
    char restart_in_filename[1024];  //!< filename for loaded restart file
    char restart_out_filename[1024]; //!< filename for saved restart file
};

/**
 * @brief All per-particle data of the simulation, stored as parallel arrays:
 * element i of each array belongs to particle i. Unlike molecular dynamics,
 * DEM particles have a radius, and they rotate: angular velocities and torques
 * are integrated alongside the translational motion.
 */
struct Vectors
{
    double time;         //!< time stamp of the vectors
    int    *type;        //!< particle types
    double *mass;        //!< masses of particles
    double *radius;      //!< radii of particles
    struct Vec3D *r;     //!< positions
    struct Vec3D *dr;    //!< displacements of the last time step
    struct Vec3D *v;     //!< velocities
    struct Vec3D *omega; //!< angular velocities
    struct Vec3D *f;     //!< forces
    struct Vec3D *T;     //!< torques
    size_t num_contacts; //!< number of particle-particle contacts, counted in @ref update_colllist
};

/**
 * @brief One entry of the neighbor list: a pair of particles within r_cut + r_shell
 * of each other, stored together with its connecting vector.
 *
 * Unlike the MD code, DEM stores the connecting vector in the pair and updates
 * it every step (@ref update_nbrlist): the contact administration reads it in
 * several places, so computing it once per step is the simpler bookkeeping.
 */
struct Pair
{
    size_t i, j;       //!< indices of the two particles forming the pair, with i < j
    struct DeltaR rij; //!< the connecting vector rij = r[i]-r[j], corrected for periodicity
};

/**
 * @brief A cell-linked-list: the box is divided into cells at least as large as
 * the interaction range, so all neighbors of a particle are found in its own
 * cell and the 26 surrounding ones.
 *
 * The particles of one cell form a linked list: head[icell] is the first
 * particle of cell icell, list[i] the particle after particle i, and SIZE_MAX
 * marks the end. This stores an arbitrary number of particles per cell in two
 * fixed-size arrays, with no per-cell allocation.
 */
struct Celllist
{
    size_t *head;                                  //!< head[icell] is the index of the first particle in cell icell, or SIZE_MAX for an empty cell
    size_t *list;                                  //!< list[i] is the next particle in the same cell as particle i; SIZE_MAX marks the end of the list
    size_t *particle2cell;                         //!< particle2cell[i] is the cell index of particle i
    size_t num_cells, num_cells_max, num_part_max; //!< number of cells used, and the numbers of cells and particles allocated for
    struct Index3D size_grid;                      //!< number of cells in each direction
};

/**
 * @brief The (Verlet) neighbor list: all particle pairs closer than
 * r_cut + r_shell when the list was last built, kept sorted by (i, j).
 *
 * The extra margin r_shell allows the same list to be reused for many steps;
 * it is rebuilt once a particle has moved more than r_shell/2 (see @ref
 * update_nbrlist). The sorting matters here: it is what lets the tangential
 * contact history tij be carried across a rebuild in a single walk over the
 * old and the new list.
 */
struct Nbrlist
{
    struct Celllist *p_celllist;   //!< cell-linked-list used to build the neighbor list
    size_t num_nbrs, num_nbrs_max; //!< number of pairs in the list, and the number allocated for
    struct Pair *nbr, *nbr_tmp;    //!< the pairs of the neighbor list, and the scratch buffer of the sorting
    /**
     * Tangential displacement of the contact between the two particles of pair
     * k, or zero while they are not touching. It lives here, beside the pair,
     * so that a contact keeps its history for as long as the pair stays in the
     * list. Only when the list is rebuilt does the history have to be carried
     * over, which @ref update_nbrlist does.
     */
    struct DeltaR *tij;
    struct Pair *nbr_prev;         //!< copy of the previous pair list, used only while the list is rebuilt
    struct DeltaR *tij_prev;       //!< copy of the previous tangential displacements, likewise
    struct DeltaR *dr;             //!< per-particle displacement since the list was last built; drives the rebuild decision
    size_t *nbr_cnt;               //!< number of pairs (i,j) per first index i, used by the counting sort in @ref build_nbrlist
};

/**
 * @brief The contacts between the particles and the walls.
 *
 * Contacts between two particles are not listed here: they are the overlapping
 * subset of the neighbor list, and the history they build up is kept beside
 * that list in @ref Nbrlist::tij. Only walls need a list of their own, since
 * they are not part of the neighbor list.
 */
struct Colllist
{
    size_t num_contacts;           //!< number of particle-particle contacts (pairs of the neighbor list that overlap)
    size_t num_contacts_max;       //!< allocated length of the contact array
    /**
     * Positions in the neighbor list of the pairs that are actually touching.
     * The overlap is tested once per time step, here, rather than again in
     * every routine that walks the contacts.
     */
    size_t *contact;
    size_t num_w;                  //!< number of particle-wall contacts
    size_t num_w_max;              //!< number of particle-wall contacts allocated for
    size_t *indcs_w;               //!< indices of the particles in contact with a wall
    unsigned int *wall_id;         //!< ID of the wall with which the particle is in contact
    struct DeltaR *riw;            //!< vector from the closest point on the wall to the particle centre, riw = ri - rw
    struct Vec3D *vw;              //!< local velocity of the wall at the contact point
    /**
     * Tangential displacement of the contact between a particle and a wall,
     * indexed as tiw[i*NUM_WALLS_MAX + w]. A particle and a wall are a fixed
     * pair, so the history of their contact can simply stay where it is instead
     * of being carried from one list of contacts to the next: it is set to zero
     * for as long as the two are not touching, and builds up while they are.
     */
    struct DeltaR *tiw;
};

#endif /* TYPES_MD_H_ */
