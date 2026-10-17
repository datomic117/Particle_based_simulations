#ifndef NBRLIST_H_
#define NBRLIST_H_

#include <stddef.h>

/**
 * @brief Allocate arrays needed to store the cell-linked-list data
 *
 * @param p_parameters used members: L, r_cut, r_shell, num_part
 * @param p_celllist the cell-linked-list to allocate
 */
void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist);

/**
 * @brief Free arrays used for the cell-linked-list
 *
 * @param p_celllist the cell-linked-list to free
 */
void free_celllist(struct Celllist *p_celllist);

/**
 * @brief Build the cell-linked-list
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_celllist 
 */
void build_celllist(struct Parameters *p_parameters, struct Vectors * p_vectors, struct Celllist *p_celllist);

/**
 * @brief Allocate arrays needed to store the neighbor list 
 * 
 * @param p_parameters 
 * @param p_nbrlist 
 */
void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist);

/**
 * @brief Free arrays use to store the neighbor list
 * 
 * @param p_nbrlist 
 */
void free_nbrlist(struct Nbrlist *p_nbrlist);

/**
 * @brief Build the neighbor list: all pairs closer than r_cut + r_shell,
 * ordered by (i, j).
 *
 * A cell-linked-list is built first, so only particles in neighboring cells
 * have to be considered as candidates; this makes the cost of the build grow
 * linearly with the number of particles instead of quadratically. The final
 * ordering allows @ref update_nbrlist to carry the tangential contact
 * histories across a rebuild in a single walk over the old and the new list.
 *
 * @param p_parameters used members: r_cut, r_shell, L
 * @param p_vectors used members: r
 * @param p_nbrlist the neighbor list to (re)build
 */
void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Comparison function used by @ref build_nbrlist to sort pairs with the
 * same first index by their second index.
 */
int cmp_sort_nbr(const void * p1, const void * p2);

/**
 * @brief Keep the neighbor list valid: rebuild it when needed, update it otherwise.
 *
 * The list contains all pairs that were closer than r_cut + r_shell when it
 * was built, and stops being trustworthy once a particle has moved more than
 * r_shell/2 since then; the list is then rebuilt, and the tangential contact
 * histories are carried over to the new list. Otherwise the stored connecting
 * vectors are updated with the displacements of the last time step.
 *
 * @param p_parameters used members: r_shell
 * @param p_vectors used members: dr
 * @param p_nbrlist used members: nbr, tij, dr, nbr_prev, tij_prev
 * @return int 1 if the list was rebuilt, 0 if it was only updated.
 */
int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Allocate memory for the collision list
 *
 * @param p_parameters used member: num_part
 * @param p_colllist the collision list to allocate
 */
void alloc_colllist(struct Parameters *p_parameters, struct Colllist *p_colllist);

/**
 * @brief Refresh the collision list: which pairs of the neighbor list overlap,
 * and which particles touch which wall. Called once per time step.
 *
 * The tangential displacement of a contact that has ended (the pair or the
 * particle and wall no longer touch) is reset to zero here.
 *
 * @param p_parameters used members: num_part, num_walls, wall_function
 * @param p_vectors used members: r, radius
 * @param p_nbrlist used members: num_nbrs, nbr, tij
 * @param p_colllist the collision list to refresh
 */
void update_colllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist,  struct Colllist* p_colllist);

/**
 * @brief Free the memory allocated for the collision list
 *
 * @param p_colllist the collision list to free
 */
void free_colllist(struct Colllist *p_colllist);

#endif /* NBRLIST */
