#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "walls.h"

void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist)
/* Allocate arrays needed to store the cell-linked-list data*/
{
    struct Index3D size_grid;
    size_t Mtot;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    const size_t num_part = p_parameters->num_part;
    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    p_celllist->head = (size_t *)malloc(Mtot * sizeof(size_t));
    p_celllist->num_cells_max = Mtot;
    p_celllist->num_cells = Mtot;
    p_celllist->num_part_max = num_part;
    p_celllist->particle2cell = (size_t *)malloc(num_part * sizeof(size_t));
    p_celllist->list = (size_t *)malloc(num_part * sizeof(size_t));
}

void free_celllist(struct Celllist *p_celllist)
/* Free arrays used for the cell-linked-list */
{
    free(p_celllist->head);
    p_celllist->head = NULL;
    free(p_celllist->particle2cell);
    p_celllist->particle2cell = NULL;
    free(p_celllist->list);
    p_celllist->list = NULL;
}

void build_celllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Celllist *p_celllist)
/* Build the cell-linked-list */
{
    struct Index3D size_grid, indx;
    size_t Mtot, icell;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    struct Vec3D mL;
    const size_t num_part = p_parameters->num_part;
    size_t *particle2cell, *head, *celllist;
    struct Vec3D *r;

    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    if (Mtot > p_celllist->num_cells_max)
    {
        p_celllist->head = (size_t *)realloc(p_celllist->head, Mtot * sizeof(size_t));
        p_celllist->num_cells_max = Mtot;
    }
    if (num_part > p_celllist->num_part_max)
    {
        p_celllist->particle2cell = (size_t *)realloc(p_celllist->particle2cell, num_part * sizeof(size_t));
        p_celllist->list = (size_t *)realloc(p_celllist->list, num_part * sizeof(size_t));
        p_celllist->num_part_max = num_part;
    }
    mL.x = ((double)size_grid.i) / p_parameters->L.x;
    mL.y = ((double)size_grid.j) / p_parameters->L.y;
    mL.z = ((double)size_grid.k) / p_parameters->L.z;
    head = p_celllist->head;
    for (icell = 0; icell < Mtot; ++icell)
        head[icell] = SIZE_MAX;
    particle2cell = p_celllist->particle2cell;
    celllist = p_celllist->list;
    r = p_vectors->r;
    // Inserting the particles in descending index order means that walking a
    // cell's linked list yields ascending indices; build_nbrlist relies on this
    // to store most pairs as (i, j) with i < j right away.
    for (size_t i = (num_part - 1); i != SIZE_MAX; --i)
    {
        indx.i = floor(r[i].x * mL.x);
        indx.j = floor(r[i].y * mL.y);
        indx.k = floor(r[i].z * mL.z);
        icell = indx.i + size_grid.i * (indx.j + indx.k * size_grid.j);
        particle2cell[i] = icell;
        celllist[i] = head[icell];
        head[icell] = i;
    }
}

void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist)
/* Allocate arrays needed to store the neighbor list */
{
    size_t num_part = p_parameters->num_part;
    double Ndouble = (double)num_part;
    double rlist = p_parameters->r_cut + p_parameters->r_shell;
    // Estimate the number of neighbors using the average density
    double Nnbr = 0.5 * Ndouble * (4.0 / 3.0 * PI * rlist * rlist * rlist) / (p_parameters->L.x * p_parameters->L.y * p_parameters->L.z);
    size_t num_nbrs_max = ceil(0.6 * Ndouble * (Nnbr + 2.0 * sqrt(Nnbr)));
    p_nbrlist->p_celllist = (struct Celllist *)malloc(sizeof(struct Celllist));
    alloc_celllist(p_parameters, p_nbrlist->p_celllist);
    p_nbrlist->num_nbrs_max = num_nbrs_max;
    p_nbrlist->num_nbrs = 0; // no pairs until build_nbrlist fills the list
    p_nbrlist->nbr = (struct Pair *)malloc(num_nbrs_max * sizeof(struct Pair));
    p_nbrlist->nbr_tmp = (struct Pair *)malloc(num_nbrs_max * sizeof(struct Pair));
    p_nbrlist->tij = (struct DeltaR *)calloc(num_nbrs_max, sizeof(struct DeltaR));
    p_nbrlist->nbr_prev = NULL;
    p_nbrlist->tij_prev = NULL;
    p_nbrlist->dr = (struct DeltaR *)malloc(p_parameters->num_part * sizeof(struct DeltaR));
    p_nbrlist->nbr_cnt = (size_t *)malloc((num_part) * sizeof(size_t));
}

void free_nbrlist(struct Nbrlist *p_nbrlist)
/* Free arrays use to store the neighbor list */
{
    free_celllist(p_nbrlist->p_celllist);
    free(p_nbrlist->p_celllist);
    p_nbrlist->p_celllist = NULL;
    free(p_nbrlist->nbr_cnt);
    p_nbrlist->nbr_cnt = NULL;
    free(p_nbrlist->nbr_tmp);
    p_nbrlist->nbr_tmp = NULL;
    free(p_nbrlist->tij);
    p_nbrlist->tij = NULL;
    free(p_nbrlist->nbr_prev);
    p_nbrlist->nbr_prev = NULL;
    free(p_nbrlist->tij_prev);
    p_nbrlist->tij_prev = NULL;
    free(p_nbrlist->nbr);
    p_nbrlist->nbr = NULL;
    free(p_nbrlist->dr);
    p_nbrlist->dr = NULL;
}

void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Build the neighbor list: all pairs closer than r_cut + r_shell, each stored
   once, ordered as (i, j) with i < j. The cell-linked-list restricts the
   candidates to particles in the same or adjacent cells, which makes the cost
   of the build grow linearly with the number of particles. The list ends up
   sorted by i and, within one i, by j — the order update_nbrlist relies on to
   carry the tangential contact history across a rebuild. */
{
    struct Index3D indx, indx_nbr;
    size_t icell, inbr;
    size_t num_nbrs, num_nbrs_max;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell; /* the radius for inclusion in the list is r_cut + r_shell */
    const double rlist_sq = rlist * rlist;
    struct Vec3D ri;
    struct Vec3D *r = p_vectors->r;
    struct DeltaR rij;
    const struct DeltaR dr = {0.0, 0.0, 0.0, 0.0};
    size_t *head, *particle2cell, *celllist;
    const int nbr_indcs[13][3] = {{0, 0, 1}, {0, 1, -1}, {0, 1, 0}, {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};
    size_t num_part = p_parameters->num_part;

    // First build a cell-linked-list
    build_celllist(p_parameters, p_vectors, p_nbrlist->p_celllist);

    /* Use the cell-linked-list to build the neighbor list: a pair is included
       when its distance is less than r_cut + r_shell. */
    struct Index3D size_grid = p_nbrlist->p_celllist->size_grid;
    num_nbrs_max = p_nbrlist->num_nbrs_max;
    num_nbrs = 0;
    struct Pair *nbr = p_nbrlist->nbr;
    size_t *nbr_cnt = p_nbrlist->nbr_cnt;
    particle2cell = p_nbrlist->p_celllist->particle2cell;
    head = p_nbrlist->p_celllist->head;
    celllist = p_nbrlist->p_celllist->list;
    for (size_t i = 0; i < num_part; ++i)
        nbr_cnt[i] = 0;
    for (size_t i = 0; i < num_part; ++i)
    {
        // Find neighbors of particle i in its own cell. Starting the walk at
        // celllist[i], the particles after i in the cell's list, stores each
        // pair once (normally with j > i, see build_celllist).
        ri = r[i];
        for (size_t j = celllist[i]; j != SIZE_MAX; j = celllist[j])
        {
            rij.x = ri.x - r[j].x;
            rij.y = ri.y - r[j].y;
            rij.z = ri.z - r[j].z;
            rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
            if (rij.sq < rlist_sq)
            {
                if (num_nbrs >= num_nbrs_max)
                {
                    num_nbrs_max += 5 * p_parameters->num_part;
                    nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                    p_nbrlist->nbr = nbr;
                    p_nbrlist->num_nbrs_max = num_nbrs_max;
                }
                if (j > i) //expected to be true unless new particles have been inserted
                {
                    nbr[num_nbrs].i = i;
                    nbr[num_nbrs].j = j;
                    ++(nbr_cnt[i]);
                }
                else //swap particles i and j
                {
                    nbr[num_nbrs].i = j;
                    nbr[num_nbrs].j = i;
                    rij.x = -rij.x;
                    rij.y = -rij.y;
                    rij.z = -rij.z;
                    ++(nbr_cnt[j]);
                }
                nbr[num_nbrs].rij = rij;
                num_nbrs++;
            }
        }

        // Next find neighbors of particle i in the 13 cells of its half-stencil:
        // visiting only half of the 26 neighboring cells stores each pair once,
        // since the pair with j in the opposite cell is found when the loop
        // reaches particle j.
        icell = particle2cell[i];
        indx.i = icell % size_grid.i;
        icell = icell / size_grid.i;
        indx.j = icell % size_grid.j;
        indx.k = icell / size_grid.j;
        for (size_t k = 0; k < 13; ++k)
        {
            ri = r[i];
            indx_nbr.i = (indx.i + nbr_indcs[k][0]);
            indx_nbr.j = (indx.j + nbr_indcs[k][1]);
            indx_nbr.k = (indx.k + nbr_indcs[k][2]);
            // The if-statements below implement the periodic boundary
            // conditions: a neighbor cell beyond the edge of the grid wraps to
            // the other side, and the position of particle i is shifted by a
            // box length so plain subtraction gives the correct distance to
            // the particles in that cell. (Note that an unsigned index that
            // would go negative wraps to SIZE_MAX.)
            if (indx_nbr.i == SIZE_MAX)
            {
                ri.x += p_parameters->L.x;
                indx_nbr.i = size_grid.i - 1;
            }
            else if (indx_nbr.i >= size_grid.i)
            {
                ri.x -= p_parameters->L.x;
                indx_nbr.i = 0;
            }
            if (indx_nbr.j == SIZE_MAX)
            {
                ri.y += p_parameters->L.y;
                indx_nbr.j = size_grid.j - 1;
            }
            else if (indx_nbr.j >= size_grid.j)
            {
                ri.y -= p_parameters->L.y;
                indx_nbr.j = 0;
            }
            if (indx_nbr.k == SIZE_MAX)
            {
                ri.z += p_parameters->L.z;
                indx_nbr.k = size_grid.k - 1;
            }
            else if (indx_nbr.k >= size_grid.k)
            {
                ri.z -= p_parameters->L.z;
                indx_nbr.k = 0;
            }
            inbr = indx_nbr.i + size_grid.i * (indx_nbr.j + indx_nbr.k * size_grid.j);
            for (size_t j = head[inbr]; j != SIZE_MAX; j = celllist[j])
            {
                rij.x = ri.x - r[j].x;
                rij.y = ri.y - r[j].y;
                rij.z = ri.z - r[j].z;
                rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
                if (rij.sq < rlist_sq)
                {
                    if (num_nbrs >= num_nbrs_max)
                    {
                        num_nbrs_max += 5 * p_parameters->num_part;
                        nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                        p_nbrlist->nbr = nbr;
                        p_nbrlist->num_nbrs_max = num_nbrs_max;
                    }
                    if (j > i)
                    {
                        nbr[num_nbrs].i = i;
                        nbr[num_nbrs].j = j;
                        ++(nbr_cnt[i]);
                    }
                    else //swap particles i and j
                    {
                        nbr[num_nbrs].i = j;
                        nbr[num_nbrs].j = i;
                        rij.x = -rij.x;
                        rij.y = -rij.y;
                        rij.z = -rij.z;
                        ++(nbr_cnt[j]);
                    }
                    nbr[num_nbrs].rij = rij;
                    num_nbrs++;
                }
            }
        }
    }

    // Sort the list by first index i with a counting sort: nbr_cnt[i] holds how
    // many pairs have first index i, the running sum below turns it into the
    // position where the pairs of each i start, and the copy that follows drops
    // every pair into its slot. A qsort per i-block then orders the second
    // index, leaving the whole list ordered by (i, j).
    size_t cnt_sum = 0;
    for (size_t i = 0; i < num_part; ++i)
    {
        size_t tmp = nbr_cnt[i];
        nbr_cnt[i] = cnt_sum;
        cnt_sum += tmp;
    }
    p_nbrlist->nbr_tmp = (struct Pair *)realloc(p_nbrlist->nbr_tmp, num_nbrs_max * sizeof(struct Pair));
    struct Pair *nbr_tmp = p_nbrlist->nbr_tmp;
    for (size_t k = 0; k < num_nbrs; ++k)
    {
        size_t i = (nbr[k].i < nbr[k].j ? nbr[k].i : nbr[k].j);
        nbr_tmp[nbr_cnt[i]++] = nbr[k];
    }
    p_nbrlist->nbr = nbr_tmp;

    // The tangential displacements are indexed by pair, so they follow the list:
    // every pair of the list just built starts without a contact. A caller that
    // wants to keep the history of the pairs that survived the rebuild copies it
    // back over this, see update_nbrlist. Note the bound is the new pair count
    // num_nbrs, not p_nbrlist->num_nbrs, which still holds the count of the
    // previous list until it is overwritten below: zeroing only that many entries
    // would leave a pair that the list gained reading the history of whichever
    // pair last occupied its index.
    p_nbrlist->tij = (struct DeltaR *)realloc(p_nbrlist->tij, p_nbrlist->num_nbrs_max * sizeof(struct DeltaR));
    for (size_t q = 0; q < num_nbrs; ++q)
        p_nbrlist->tij[q] = (struct DeltaR){0.0, 0.0, 0.0, 0.0};
    p_nbrlist->nbr_tmp = nbr;
    nbr = p_nbrlist->nbr;
    nbr_tmp = p_nbrlist->nbr_tmp;
    size_t k = 0;
    for (size_t i = 0; i < num_part; k = nbr_cnt[i], i++)
        qsort(nbr + k, nbr_cnt[i] - k, sizeof(struct Pair), cmp_sort_nbr);
    p_nbrlist->num_nbrs = num_nbrs;
    p_nbrlist->dr = (struct DeltaR *)realloc(p_nbrlist->dr, num_part * sizeof(struct DeltaR));
    for (size_t i = 0; i < num_part; ++i) /*initialize particle displacements (with respect to creation time) to zero */
        p_nbrlist->dr[i] = dr;
}

// Comparison function for qsort: orders pairs with the same first index by
// their second index.
int cmp_sort_nbr(const void *p1, const void *p2)
{
    const struct Pair *p_nbr1 = p1;
    const struct Pair *p_nbr2 = p2;
    if (p_nbr1->j < p_nbr2->j)
        return -1;
    else
        return (p_nbr1->j - p_nbr2->j);
}

int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Update the connecting vectors of the neighbor list and, if needed, rebuild it.
   The list contains all pairs that were within r_cut + r_shell at build time,
   so it misses a pair within r_cut only after two particles have together
   bridged the margin r_shell — which cannot happen before one of them has
   moved r_shell/2. Once a particle has, the list is rebuilt. */
{
    const double dr_sq_max = 0.25 * (p_parameters->r_shell * p_parameters->r_shell);
    int isRebuild = 0;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    struct Vec3D *dr = p_vectors->dr;
    // Check whether any particle has moved more than r_shell/2 since the build
    for (size_t i = 0; i < p_parameters->num_part; ++i)
        if ((p_nbrlist->dr[i].sq) > dr_sq_max)
        {
            isRebuild = 1;
            break; // one particle is enough, the whole list is rebuilt anyway
        }
    if (isRebuild) // rebuild neighbor list
    {
        // The pairs are about to be renumbered, so the tangential displacements
        // have to move with them. Copy the old pairs and their history aside,
        // build the new list, and carry over what the two have in common.
        // build_nbrlist finishes with a counting sort that leaves the pairs
        // ordered by i and, within one i, by j, and it swaps its own buffers
        // while doing so, which is why the copy is made rather than the buffers
        // borrowed. Both lists being ordered the same way, one walk over the two
        // finds every pair they share.
        size_t num_old = p_nbrlist->num_nbrs;
        p_nbrlist->nbr_prev = (struct Pair *)realloc(p_nbrlist->nbr_prev, num_old * sizeof(struct Pair));
        p_nbrlist->tij_prev = (struct DeltaR *)realloc(p_nbrlist->tij_prev, num_old * sizeof(struct DeltaR));
        struct Pair *nbr_old = p_nbrlist->nbr_prev;
        struct DeltaR *tij_old = p_nbrlist->tij_prev;
        for (size_t k = 0; k < num_old; ++k)
        {
            nbr_old[k] = p_nbrlist->nbr[k];
            tij_old[k] = p_nbrlist->tij[k];
        }

        build_nbrlist(p_parameters, p_vectors, p_nbrlist);

        size_t num_new = p_nbrlist->num_nbrs;
        struct Pair *nbr_new = p_nbrlist->nbr;
        struct DeltaR *tij_new = p_nbrlist->tij;

        for (size_t m = 0, k = 0; m < num_new && k < num_old;)
        {
            if (nbr_new[m].i != nbr_old[k].i)
                (nbr_new[m].i < nbr_old[k].i) ? ++m : ++k;
            else if (nbr_new[m].j != nbr_old[k].j)
                (nbr_new[m].j < nbr_old[k].j) ? ++m : ++k;
            else // the same pair is in both lists, so its contact continues
                tij_new[m++] = tij_old[k++];
        }
    }
    else // If no rebuild is needed, update the stored connecting vectors by the
         // displacements both particles made this time step
    {
        for (size_t k = 0; k < p_nbrlist->num_nbrs; ++k)
        {
            size_t i = nbr[k].i;
            size_t j = nbr[k].j;
            rij = nbr[k].rij;
            rij.x += (dr[i].x - dr[j].x);
            rij.y += (dr[i].y - dr[j].y);
            rij.z += (dr[i].z - dr[j].z);
            rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
            p_nbrlist->nbr[k].rij = rij;
        }
    }
    return isRebuild;
}

void alloc_colllist(struct Parameters *p_parameters, struct Colllist *p_colllist)
{
    p_colllist->num_contacts = 0;
    p_colllist->num_contacts_max = 0;
    p_colllist->contact = NULL;
    p_colllist->num_w = 0;
    size_t num_w_max = 0;
    p_colllist->num_w_max = num_w_max;
    p_colllist->indcs_w = (size_t *)malloc(num_w_max * sizeof(size_t));
    p_colllist->wall_id = (unsigned int *)malloc(num_w_max * sizeof(unsigned int));
    p_colllist->riw = (struct DeltaR *)malloc(num_w_max * sizeof(struct DeltaR));
    p_colllist->vw = (struct Vec3D *)malloc(num_w_max * sizeof(struct Vec3D));
    // One tangential displacement per particle and wall, held for the whole run
    // and zero whenever the two are not in contact.
    p_colllist->tiw = (struct DeltaR *)calloc(p_parameters->num_part * NUM_WALLS_MAX, sizeof(struct DeltaR));
}

void update_colllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist)
/* Collects the contacts between the particles and the walls.
   Contacts between two particles are not collected: they are the overlapping
   subset of the neighbor list, which calculate_forces_pp and
   update_tangential_displacements walk directly. */
{
    double * R = p_vectors->radius;

    // Collect the pairs of the neighbor list that overlap. A pair that is not
    // touching has no contact, so whatever tangential displacement it had is
    // discarded here.
    size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Pair *nbr = p_nbrlist->nbr;
    struct DeltaR *tij = p_nbrlist->tij;
    if (num_nbrs > p_colllist->num_contacts_max)
    {
        p_colllist->num_contacts_max = num_nbrs;
        p_colllist->contact = (size_t *)realloc(p_colllist->contact, num_nbrs * sizeof(size_t));
    }
    size_t *contact = p_colllist->contact;
    size_t num_contacts = 0;
    for (size_t k = 0; k < num_nbrs; ++k)
    {
        double sumR = R[nbr[k].i] + R[nbr[k].j];
        if (nbr[k].rij.sq < sumR * sumR)
            contact[num_contacts++] = k;
        else
            tij[k] = (struct DeltaR){0.0, 0.0, 0.0, 0.0};
    }
    p_colllist->num_contacts = num_contacts;

    /* Determine the particles in collision with a wall */
    size_t num_part = p_parameters->num_part;
    size_t *indcs_w = p_colllist->indcs_w;
    unsigned int *wall_id = p_colllist->wall_id;
    struct DeltaR *tiw = p_colllist->tiw;
    struct DeltaR *riw = p_colllist->riw;
    struct Vec3D *vw = p_colllist->vw;
    struct Vec3D *r = p_vectors->r;
    unsigned int num_walls = p_parameters->num_walls;
    size_t num_w_max = p_colllist->num_w_max;
    size_t k = 0;
    for (size_t i = 0; i < num_part; i++)
    {
        for (unsigned int j = 0; j < num_walls; ++j)
        {
            struct DeltaR riw_loc;
            struct Vec3D vw_loc;
            if (!p_parameters->wall_function[j](p_parameters, R[i], &r[i], &riw_loc, &vw_loc))
            {
                // Not touching: the contact, if there was one, is over and its
                // tangential displacement goes back to zero.
                tiw[i * NUM_WALLS_MAX + j] = (struct DeltaR){0.0, 0.0, 0.0, 0.0};
            }
            else
            {
                if (k >= num_w_max)
                {
                    num_w_max += 10;
                    indcs_w = (size_t *)realloc(indcs_w, num_w_max * sizeof(size_t));
                    wall_id = (unsigned int *)realloc(wall_id, num_w_max * sizeof(unsigned int));
                    riw = (struct DeltaR *)realloc(riw, num_w_max * sizeof(struct DeltaR));
                    vw = (struct Vec3D *)realloc(vw, num_w_max * sizeof(struct Vec3D));
                }
                indcs_w[k] = i;
                wall_id[k] = j;
                riw[k] = riw_loc;
                vw[k] = vw_loc;
                ++k;
            }
        }
    }
    size_t num_w = k;
    p_colllist->num_w = num_w;
    p_colllist->num_w_max = num_w_max;
    p_colllist->indcs_w = indcs_w;
    p_colllist->wall_id = wall_id;
    p_colllist->riw = riw;
    p_colllist->vw = vw;
}

void free_colllist(struct Colllist *p_colllist)
{
    free(p_colllist->contact);
    free(p_colllist->indcs_w);
    free(p_colllist->wall_id);
    free(p_colllist->riw);
    free(p_colllist->tiw);
    free(p_colllist->vw);
}
