#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include "pairlist.h"

#define ADDRESS(x,y,z) (((z)*GY + (y))*GX + (x))

static int pairlist1(int nAtoms, PL_FLOAT *atoms, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);
static int pairlist2(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);

/* Return 0 if ngrid[d] >= 1 for all d, else -1. */
static int check_ngrid(const int ngrid[3])
{
  for (int d = 0; d < 3; d++) {
    if (ngrid[d] < 1) {
      return -1;
    }
  }
  return 0;
}

/* GX*GY*GZ must fit in int (ADDRESS and indexing use int). */
static int compute_n_total_grids(const int ngrid[3], int *nTotalGrids)
{
  long long t = (long long)ngrid[0] * (long long)ngrid[1];
  if (t > INT_MAX) {
    return -1;
  }
  t *= (long long)ngrid[2];
  if (t > INT_MAX) {
    return -1;
  }
  *nTotalGrids = (int)t;
  return 0;
}

/* *sum += a * b; return -1 on signed int overflow. */
static int accum_prod(int *sum, int a, int b)
{
  if (a < 0 || b < 0) {
    return -1;
  }
  if (a > 0 && b > INT_MAX / a) {
    return -1;
  }
  int prod = a * b;
  if (*sum > INT_MAX - prod) {
    return -1;
  }
  *sum += prod;
  return 0;
}

/* Allocate n ints; n==0 → NULL (empty). Overflow/OOM → NULL. */
static int *malloc_int_array(size_t n)
{
  if (n == 0) {
    return NULL;
  }
  if (n > SIZE_MAX / sizeof(int)) {
    return NULL;
  }
  return (int *)malloc(n * sizeof(int));
}

/* Map a fractional coordinate to a grid index in [0, ngrid_d - 1].
   For finite x, x-=floor(x) yields [0,1); under IEEE-754, floor(x*ngrid_d)
   is then always in [0, ngrid_d-1], so no clamp is required. */
static inline int frac_to_grid(PL_FLOAT x, int ngrid_d)
{
  x -= floor(x);
  return (int)floor(x * (PL_FLOAT)ngrid_d);
}

/* Neighbor offset range along one axis of size G.
   G==1: {0}; G==2: {0,+1} (avoids -1≡+1 alias); G>=3: {-1,0,+1}. */
static inline void neighbor_range(int G, int i, int *k0, int *k1)
{
  *k0 = (G < 3) ? i : i - 1;
  *k1 = (G < 2) ? i : i + 1;
}

typedef struct {
  int *nResidents;
  int *residents;
  int *heads;
} ResidentList;

static void resident_list_clear(ResidentList *rl)
{
  free(rl->nResidents);
  free(rl->residents);
  free(rl->heads);
  rl->nResidents = NULL;
  rl->residents = NULL;
  rl->heads = NULL;
}

/* Bin atoms into grid cells. On failure leaves *rl cleared and returns -1. */
static int resident_list_build(int npos, const PL_FLOAT *rpos, const int ngrid[3],
                               int nTotalGrids, ResidentList *rl)
{
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  int *pointer = NULL;

  rl->nResidents = NULL;
  rl->residents = NULL;
  rl->heads = NULL;

  rl->nResidents = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (rl->nResidents == NULL) {
    return -1;
  }
  for (int i = 0; i < npos; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos[i * 3 + d], ngrid[d]);
    }
    rl->nResidents[ADDRESS(grid[0], grid[1], grid[2])]++;
  }

  if ((size_t)npos > SIZE_MAX - (size_t)nTotalGrids) {
    resident_list_clear(rl);
    return -1;
  }
  pointer = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  rl->residents = (int *)calloc((size_t)npos + (size_t)nTotalGrids, sizeof(int));
  rl->heads = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (pointer == NULL || rl->residents == NULL || rl->heads == NULL) {
    free(pointer);
    resident_list_clear(rl);
    return -1;
  }

  int head = 0;
  for (int g = 0; g < nTotalGrids; g++) {
    rl->heads[g] = head;
    pointer[g] = head;
    head += rl->nResidents[g] + 1;
    rl->residents[head - 1] = -1;
  }
  for (int i = 0; i < npos; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos[i * 3 + d], ngrid[d]);
    }
    int a = ADDRESS(grid[0], grid[1], grid[2]);
    rl->residents[pointer[a]] = i;
    pointer[a]++;
  }
  free(pointer);
  return 0;
}

/* Shrink malloc'd pair buffer from estim pairs to nPairs. Keeps buffer if realloc fails. */
static void shrink_pair_buffer(int **pairs, int estim, int nPairs)
{
  if (nPairs == 0) {
    free(*pairs);
    *pairs = NULL;
    return;
  }
  if (nPairs >= estim) {
    return;
  }
  int *p = (int *)realloc(*pairs, (size_t)nPairs * 2u * sizeof(int));
  if (p != NULL) {
    *pairs = p;
  }
}

int
pairlist(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs)
{
  if ((nAtoms1 == 0) || (atoms1 == NULL) || (atoms0 == atoms1)) {
    return pairlist1(nAtoms0, atoms0, lower, higher, cell, pairs);
  }
  return pairlist2(nAtoms0, atoms0, nAtoms1, atoms1, lower, higher, cell, pairs);
}



static int
gridpairlist(int ngrid[3], int single, int **pairs)
/* 
given:
    ngrid: number of grids in three dimensions
    single: if 1, make pair list for a single component system
returns:
    return value: number of grid pairs
    pairs: newly allocated list of grid pairs
*/
{
  if (check_ngrid(ngrid) < 0) {
    *pairs = NULL;
    return -1;
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    *pairs = NULL;
    return -1;
  }
  int npair = 0;
  int neighborcells = 27;
  if ( single ) neighborcells = 14;
  /* neighborcells * nTotalGrids * 2 ints */
  if (nTotalGrids > INT_MAX / (neighborcells * 2)) {
    *pairs = NULL;
    return -1;
  }
  (*pairs) = malloc_int_array((size_t)neighborcells * (size_t)nTotalGrids * 2u);
  if ((*pairs) == NULL) {
    return -1;
  }
  //make neighboring grid list
  for(int ix=0;ix<GX;ix++){
    for(int iy=0;iy<GY;iy++){
      for(int iz=0;iz<GZ;iz++){
        int a = ADDRESS(ix,iy,iz); //center
        int kx0, kx1, ky0, ky1, kz0, kz1;
        neighbor_range(GX, ix, &kx0, &kx1);
        for(int kx=kx0;kx<=kx1;kx++){
          int jx = (kx+GX) % GX;
          neighbor_range(GY, iy, &ky0, &ky1);
          for(int ky=ky0;ky<=ky1;ky++){
            int jy = (ky+GY) % GY;
            neighbor_range(GZ, iz, &kz0, &kz1);
            for(int kz=kz0;kz<=kz1;kz++){
              int jz = (kz+GZ) % GZ;
              int b = ADDRESS(jx,jy,jz);
              if ( ( ! single ) || ( a<b ) ){
                (*pairs)[npair*2+0] = a;
                (*pairs)[npair*2+1] = b;
                npair ++;
              }
            }
          }
        }
      }
    }
  }
  return npair;
}




//Returns the rough neighbor list
int
Pairs(int npos, PL_FLOAT* rpos, int ngrid[3], int **pairs)
{
/*
given:
    npos: number of atoms
    rpos: fractional positions of atoms
    ngrid: number of grid divisions
returns:
    return value: number of pairs (>=0), or -1 on error
    pairs: newly allocated list of pairs (NULL if count is 0 or on error)
*/
  ResidentList rl = {0};
  int *gridPairs = NULL;
  *pairs = NULL;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }
  if (resident_list_build(npos, rpos, ngrid, nTotalGrids, &rl) < 0) {
    return -1;
  }

  int nGridPairs = gridpairlist(ngrid, 1, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  int estim = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    if (accum_prod(&estim, rl.nResidents[g0], rl.nResidents[g1]) < 0) {
      goto fail;
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    int n = rl.nResidents[g];
    if (n > 1) {
      if (n - 1 > INT_MAX / n) {
        goto fail;
      }
      int half = n * (n - 1) / 2;
      if (estim > INT_MAX - half) {
        goto fail;
      }
      estim += half;
    }
  }

  if (estim > INT_MAX / 2) {
    goto fail;
  }
  *pairs = malloc_int_array((size_t)estim * 2u);
  if (estim > 0 && *pairs == NULL) {
    goto fail;
  }

  int nPairs = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl.nResidents[g0]; j++) {
      int r0 = rl.residents[rl.heads[g0] + j];
      for (int k = 0; k < rl.nResidents[g1]; k++) {
        int r1 = rl.residents[rl.heads[g1] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    for (int j = 0; j < rl.nResidents[g]; j++) {
      int r0 = rl.residents[rl.heads[g] + j];
      for (int k = j + 1; k < rl.nResidents[g]; k++) {
        int r1 = rl.residents[rl.heads[g] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }

  shrink_pair_buffer(pairs, estim, nPairs);
  free(gridPairs);
  resident_list_clear(&rl);
  return nPairs;

fail:
  free(gridPairs);
  resident_list_clear(&rl);
  free(*pairs);
  *pairs = NULL;
  return -1;
}



int
Pairs2(int npos0, PL_FLOAT *rpos0, int npos1, PL_FLOAT *rpos1, int ngrid[3], int **pairs)
/*
given:
    npos0: number of atoms of component 0
    rpos0: fractional positions of atoms
    npos1: number of atoms of component 1
    rpos1: fractional positions of atoms
    ngrid: grid divisions of the cell
returns:
    return value: number of pairs (>=0), or -1 on error
    pairs: newly allocated list of pairs (NULL if count is 0 or on error)
*/
{
  ResidentList rl0 = {0};
  ResidentList rl1 = {0};
  int *gridPairs = NULL;
  *pairs = NULL;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }
  if (resident_list_build(npos0, rpos0, ngrid, nTotalGrids, &rl0) < 0) {
    return -1;
  }
  if (resident_list_build(npos1, rpos1, ngrid, nTotalGrids, &rl1) < 0) {
    resident_list_clear(&rl0);
    return -1;
  }

  int nGridPairs = gridpairlist(ngrid, 0, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  /* gridPairs already includes same-cell (g,g); do not count them again. */
  int estim = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    if (accum_prod(&estim, rl0.nResidents[g0], rl1.nResidents[g1]) < 0) {
      goto fail;
    }
  }

  if (estim > INT_MAX / 2) {
    goto fail;
  }
  *pairs = malloc_int_array((size_t)estim * 2u);
  if (estim > 0 && *pairs == NULL) {
    goto fail;
  }

  int nPairs = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl0.nResidents[g0]; j++) {
      int r0 = rl0.residents[rl0.heads[g0] + j];
      for (int k = 0; k < rl1.nResidents[g1]; k++) {
        int r1 = rl1.residents[rl1.heads[g1] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }

  shrink_pair_buffer(pairs, estim, nPairs);
  free(gridPairs);
  resident_list_clear(&rl0);
  resident_list_clear(&rl1);
  return nPairs;

fail:
  free(gridPairs);
  resident_list_clear(&rl0);
  resident_list_clear(&rl1);
  free(*pairs);
  *pairs = NULL;
  return -1;
}



/* Minimum-image squared distance for fractional coords + triclinic cell
   (rows a,b,c in cell[9] row-major). Matches np.dot(d, cell) then ||.||^2. */
static inline PL_FLOAT frac_minimage_dist2(const PL_FLOAT *ri, const PL_FLOAT *rj,
                                           const PL_FLOAT cell[9])
{
  PL_FLOAT dx = ri[0] - rj[0];
  PL_FLOAT dy = ri[1] - rj[1];
  PL_FLOAT dz = ri[2] - rj[2];
  dx -= floor(dx + 0.5);
  dy -= floor(dy + 0.5);
  dz -= floor(dz + 0.5);
  PL_FLOAT rx = dx * cell[0] + dy * cell[3] + dz * cell[6];
  PL_FLOAT ry = dx * cell[1] + dy * cell[4] + dz * cell[7];
  PL_FLOAT rz = dx * cell[2] + dy * cell[5] + dz * cell[8];
  return rx * rx + ry * ry + rz * rz;
}


int
PairsFiltered(int npos, PL_FLOAT *rpos, int ngrid[3], const PL_FLOAT cell[9],
              PL_FLOAT rc, int **pairs, PL_FLOAT **dists)
{
  ResidentList rl = {0};
  int *gridPairs = NULL;
  *pairs = NULL;
  if (dists) {
    *dists = NULL;
  }
  if (rc < 0.0) {
    return -1;
  }
  const PL_FLOAT rc2 = rc * rc;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }
  if (resident_list_build(npos, rpos, ngrid, nTotalGrids, &rl) < 0) {
    return -1;
  }

  int nGridPairs = gridpairlist(ngrid, 1, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  /* Pass 1: count pairs with dist < rc */
  int nPairs = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl.nResidents[g0]; j++) {
      int r0 = rl.residents[rl.heads[g0] + j];
      for (int k = 0; k < rl.nResidents[g1]; k++) {
        int r1 = rl.residents[rl.heads[g1] + k];
        if (frac_minimage_dist2(rpos + r0 * 3, rpos + r1 * 3, cell) < rc2) {
          nPairs++;
        }
      }
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    for (int j = 0; j < rl.nResidents[g]; j++) {
      int r0 = rl.residents[rl.heads[g] + j];
      for (int k = j + 1; k < rl.nResidents[g]; k++) {
        int r1 = rl.residents[rl.heads[g] + k];
        if (frac_minimage_dist2(rpos + r0 * 3, rpos + r1 * 3, cell) < rc2) {
          nPairs++;
        }
      }
    }
  }

  *pairs = malloc_int_array((size_t)nPairs * 2u);
  if (nPairs > 0 && *pairs == NULL) {
    goto fail;
  }
  if (dists) {
    if (nPairs == 0) {
      *dists = NULL;
    } else {
      *dists = (PL_FLOAT *)malloc(sizeof(PL_FLOAT) * (size_t)nPairs);
      if (*dists == NULL) {
        goto fail;
      }
    }
  }

  /* Pass 2: write */
  int store = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl.nResidents[g0]; j++) {
      int r0 = rl.residents[rl.heads[g0] + j];
      for (int k = 0; k < rl.nResidents[g1]; k++) {
        int r1 = rl.residents[rl.heads[g1] + k];
        PL_FLOAT d2 = frac_minimage_dist2(rpos + r0 * 3, rpos + r1 * 3, cell);
        if (d2 < rc2) {
          (*pairs)[store * 2 + 0] = r0;
          (*pairs)[store * 2 + 1] = r1;
          if (dists) {
            (*dists)[store] = sqrt(d2);
          }
          store++;
        }
      }
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    for (int j = 0; j < rl.nResidents[g]; j++) {
      int r0 = rl.residents[rl.heads[g] + j];
      for (int k = j + 1; k < rl.nResidents[g]; k++) {
        int r1 = rl.residents[rl.heads[g] + k];
        PL_FLOAT d2 = frac_minimage_dist2(rpos + r0 * 3, rpos + r1 * 3, cell);
        if (d2 < rc2) {
          (*pairs)[store * 2 + 0] = r0;
          (*pairs)[store * 2 + 1] = r1;
          if (dists) {
            (*dists)[store] = sqrt(d2);
          }
          store++;
        }
      }
    }
  }

  free(gridPairs);
  resident_list_clear(&rl);
  return store;

fail:
  free(gridPairs);
  resident_list_clear(&rl);
  free(*pairs);
  *pairs = NULL;
  if (dists) {
    free(*dists);
    *dists = NULL;
  }
  return -1;
}


int
Pairs2Filtered(int npos0, PL_FLOAT *rpos0, int npos1, PL_FLOAT *rpos1,
               int ngrid[3], const PL_FLOAT cell[9], PL_FLOAT rc,
               int **pairs, PL_FLOAT **dists)
{
  ResidentList rl0 = {0};
  ResidentList rl1 = {0};
  int *gridPairs = NULL;
  *pairs = NULL;
  if (dists) {
    *dists = NULL;
  }
  if (rc < 0.0) {
    return -1;
  }
  const PL_FLOAT rc2 = rc * rc;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }
  if (resident_list_build(npos0, rpos0, ngrid, nTotalGrids, &rl0) < 0) {
    return -1;
  }
  if (resident_list_build(npos1, rpos1, ngrid, nTotalGrids, &rl1) < 0) {
    resident_list_clear(&rl0);
    return -1;
  }

  int nGridPairs = gridpairlist(ngrid, 0, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  int nPairs = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl0.nResidents[g0]; j++) {
      int r0 = rl0.residents[rl0.heads[g0] + j];
      for (int k = 0; k < rl1.nResidents[g1]; k++) {
        int r1 = rl1.residents[rl1.heads[g1] + k];
        if (frac_minimage_dist2(rpos0 + r0 * 3, rpos1 + r1 * 3, cell) < rc2) {
          nPairs++;
        }
      }
    }
  }

  *pairs = malloc_int_array((size_t)nPairs * 2u);
  if (nPairs > 0 && *pairs == NULL) {
    goto fail;
  }
  if (dists) {
    if (nPairs == 0) {
      *dists = NULL;
    } else {
      *dists = (PL_FLOAT *)malloc(sizeof(PL_FLOAT) * (size_t)nPairs);
      if (*dists == NULL) {
        goto fail;
      }
    }
  }

  int store = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    for (int j = 0; j < rl0.nResidents[g0]; j++) {
      int r0 = rl0.residents[rl0.heads[g0] + j];
      for (int k = 0; k < rl1.nResidents[g1]; k++) {
        int r1 = rl1.residents[rl1.heads[g1] + k];
        PL_FLOAT d2 = frac_minimage_dist2(rpos0 + r0 * 3, rpos1 + r1 * 3, cell);
        if (d2 < rc2) {
          (*pairs)[store * 2 + 0] = r0;
          (*pairs)[store * 2 + 1] = r1;
          if (dists) {
            (*dists)[store] = sqrt(d2);
          }
          store++;
        }
      }
    }
  }

  free(gridPairs);
  resident_list_clear(&rl0);
  resident_list_clear(&rl1);
  return store;

fail:
  free(gridPairs);
  resident_list_clear(&rl0);
  resident_list_clear(&rl1);
  free(*pairs);
  *pairs = NULL;
  if (dists) {
    free(*dists);
    *dists = NULL;
  }
  return -1;
}




static int
pairlist1(int nAtoms, PL_FLOAT *atoms, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs)
/* 
given:
    nAtoms: number of atoms
    atoms: positions of atoms
    lower, higher: two atoms are considered to be paired when their distance is between them.
    cell: rectangular cell shape
returns:
    return value: number of pairs
    pairs: newly allocated list of pairs
*/
{
  *pairs = NULL;
  if (nAtoms > 0 && (size_t)nAtoms > SIZE_MAX / (sizeof(PL_FLOAT) * 3u)) {
    return -1;
  }
  PL_FLOAT *rpos = NULL;
  if (nAtoms > 0) {
    rpos = (PL_FLOAT *)malloc(sizeof(PL_FLOAT) * 3u * (size_t)nAtoms);
    if (rpos == NULL) {
      return -1;
    }
  }
  for(int i=0; i<nAtoms; i++){
    for(int d=0; d<3; d++){
      PL_FLOAT r = atoms[i*3+d] / cell[d];
      r -= floor(r);
      rpos[i*3+d] = r;
    }
  }
  int ngrid[3];
  for(int d=0;d<3;d++){
    ngrid[d] = (int) floor(cell[d] / higher);
    if (ngrid[d] < 1) {
      ngrid[d] = 1;
    }
  }
  int npairs = Pairs(nAtoms, rpos, ngrid, pairs);
  if (npairs < 0) {
    free(rpos);
    return -1;
  }
  int store=0;
  for(int head=0; head<npairs; head++){
    int r0 = (*pairs)[head*2];
    int r1 = (*pairs)[head*2+1];
    PL_FLOAT sum = 0.0;
    for(int d=0;d<3;d++){
      PL_FLOAT delta = rpos[r0*3+d] - rpos[r1*3+d];
      delta -= floor(delta+0.5);
      delta *= cell[d];
      sum += delta*delta;
    }
    if ((lower*lower < sum) && (sum < higher*higher)){
      (*pairs)[store*2+0] = r0;
      (*pairs)[store*2+1] = r1;
      store ++;
    }
  }
  free(rpos);
  return store;
}



static int
pairlist2(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs)
/* 
given:
    nAtoms0: number of atoms of component 0
    atoms0: positions of atoms
    nAtoms1: number of atoms of component 1
    atoms1: positions of atoms
    higher: two atoms are considered to be paired when their distance is below it.
    cell: rectangular cell shape
returns:
    return value: number of pairs
    pairs: newly allocated list of pairs
*/
{
  *pairs = NULL;
  if (nAtoms0 > 0 && (size_t)nAtoms0 > SIZE_MAX / (sizeof(PL_FLOAT) * 3u)) {
    return -1;
  }
  if (nAtoms1 > 0 && (size_t)nAtoms1 > SIZE_MAX / (sizeof(PL_FLOAT) * 3u)) {
    return -1;
  }
  PL_FLOAT *rpos0 = NULL;
  PL_FLOAT *rpos1 = NULL;
  if (nAtoms0 > 0) {
    rpos0 = (PL_FLOAT *)malloc(sizeof(PL_FLOAT) * 3u * (size_t)nAtoms0);
    if (rpos0 == NULL) {
      return -1;
    }
  }
  for(int i=0; i<nAtoms0; i++){
    for(int d=0; d<3; d++){
      PL_FLOAT r = atoms0[i*3+d] / cell[d];
      r -= floor(r);
      rpos0[i*3+d] = r;
    }
  }
  if (nAtoms1 > 0) {
    rpos1 = (PL_FLOAT *)malloc(sizeof(PL_FLOAT) * 3u * (size_t)nAtoms1);
    if (rpos1 == NULL) {
      free(rpos0);
      return -1;
    }
  }
  for(int i=0; i<nAtoms1; i++){
    for(int d=0; d<3; d++){
      PL_FLOAT r = atoms1[i*3+d] / cell[d];
      r -= floor(r);
      rpos1[i*3+d] = r;
    }
  }
  int ngrid[3];
  for(int d=0;d<3;d++){
    ngrid[d] = (int) floor(cell[d] / higher);
    if (ngrid[d] < 1) {
      ngrid[d] = 1;
    }
  }
  int npairs = Pairs2(nAtoms0, rpos0, nAtoms1, rpos1, ngrid, pairs);
  if (npairs < 0) {
    free(rpos0);
    free(rpos1);
    return -1;
  }
  int store=0;
  for(int head=0; head<npairs; head++){
    int r0 = (*pairs)[head*2];
    int r1 = (*pairs)[head*2+1];
    PL_FLOAT sum = 0.0;
    for(int d=0;d<3;d++){
      PL_FLOAT delta = rpos0[r0*3+d] - rpos1[r1*3+d];
      delta -= floor(delta+0.5);
      delta *= cell[d];
      sum += delta*delta;
    }
    if ((lower*lower < sum) && (sum < higher*higher)){
      (*pairs)[store*2+0] = r0;
      (*pairs)[store*2+1] = r1;
      store ++;
    }
  }
  free(rpos0);
  free(rpos1);
  return store;
}

