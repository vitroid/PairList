#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include "pairlist.h"

// neighborlist.cを使って書きなおそうかとも思ったが、効率が若干悪くなる気がするのでやめる。

#define ADDRESS(x,y,z) (((z)*GY + (y))*GX + (x))
#define True 1
#define False 0
int pairlist1(int nAtoms, PL_FLOAT *atoms, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);
int pairlist2(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);

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

int
pairlist(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs)
{
  if ( ( nAtoms1 == 0 ) || (atoms1 == NULL) || ( atoms0 == atoms1) ){
    if ( atoms0 == atoms1)
      assert (nAtoms0 == nAtoms1);
    return pairlist1(nAtoms0, atoms0, lower, higher, cell, pairs);
  }
  return pairlist2(nAtoms0, atoms0, nAtoms1, atoms1, lower, higher, cell, pairs);
}
   


int
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
  int *nResidents = NULL;
  int *pointer = NULL;
  int *residents = NULL;
  int *heads = NULL;
  int *gridPairs = NULL;
  *pairs = NULL;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }

  nResidents = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (nResidents == NULL) {
    goto fail;
  }
  for (int i = 0; i < npos; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos[i * 3 + d], ngrid[d]);
    }
    nResidents[ADDRESS(grid[0], grid[1], grid[2])]++;
  }

  if ((size_t)npos > SIZE_MAX - (size_t)nTotalGrids) {
    goto fail;
  }
  pointer = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  residents = (int *)calloc((size_t)npos + (size_t)nTotalGrids, sizeof(int));
  heads = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (pointer == NULL || residents == NULL || heads == NULL) {
    goto fail;
  }
  int head = 0;
  for (int g = 0; g < nTotalGrids; g++) {
    heads[g] = head;
    pointer[g] = head;
    head += nResidents[g] + 1;
    residents[head - 1] = -1; /* terminator */
  }
  for (int i = 0; i < npos; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos[i * 3 + d], ngrid[d]);
    }
    int a = ADDRESS(grid[0], grid[1], grid[2]);
    residents[pointer[a]] = i;
    pointer[a]++;
  }

  int nGridPairs = gridpairlist(ngrid, True, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  int estim = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    if (accum_prod(&estim, nResidents[g0], nResidents[g1]) < 0) {
      goto fail;
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    int n = nResidents[g];
    if (n > 1) {
      /* n*(n-1)/2 without overflowing the addend */
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
    for (int j = 0; j < nResidents[g0]; j++) {
      int r0 = residents[heads[g0] + j];
      for (int k = 0; k < nResidents[g1]; k++) {
        int r1 = residents[heads[g1] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }
  for (int g = 0; g < nTotalGrids; g++) {
    for (int j = 0; j < nResidents[g]; j++) {
      int r0 = residents[heads[g] + j];
      for (int k = j + 1; k < nResidents[g]; k++) {
        int r1 = residents[heads[g] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }

  free(gridPairs);
  free(nResidents);
  free(pointer);
  free(residents);
  free(heads);
  return nPairs;

fail:
  free(gridPairs);
  free(nResidents);
  free(pointer);
  free(residents);
  free(heads);
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
  int *nResidents0 = NULL;
  int *nResidents1 = NULL;
  int *pointer = NULL;
  int *residents0 = NULL;
  int *heads0 = NULL;
  int *residents1 = NULL;
  int *heads1 = NULL;
  int *gridPairs = NULL;
  *pairs = NULL;

  if (check_ngrid(ngrid) < 0) {
    return -1;
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  int nTotalGrids;
  if (compute_n_total_grids(ngrid, &nTotalGrids) < 0) {
    return -1;
  }

  nResidents0 = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  nResidents1 = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (nResidents0 == NULL || nResidents1 == NULL) {
    goto fail;
  }
  for (int i = 0; i < npos0; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos0[i * 3 + d], ngrid[d]);
    }
    nResidents0[ADDRESS(grid[0], grid[1], grid[2])]++;
  }
  for (int i = 0; i < npos1; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos1[i * 3 + d], ngrid[d]);
    }
    nResidents1[ADDRESS(grid[0], grid[1], grid[2])]++;
  }

  if ((size_t)npos0 > SIZE_MAX - (size_t)nTotalGrids ||
      (size_t)npos1 > SIZE_MAX - (size_t)nTotalGrids) {
    goto fail;
  }
  pointer = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  residents0 = (int *)calloc((size_t)npos0 + (size_t)nTotalGrids, sizeof(int));
  heads0 = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (pointer == NULL || residents0 == NULL || heads0 == NULL) {
    goto fail;
  }
  int head = 0;
  for (int g = 0; g < nTotalGrids; g++) {
    heads0[g] = head;
    pointer[g] = head;
    head += nResidents0[g] + 1;
    residents0[head - 1] = -1;
  }
  for (int i = 0; i < npos0; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos0[i * 3 + d], ngrid[d]);
    }
    int a = ADDRESS(grid[0], grid[1], grid[2]);
    residents0[pointer[a]] = i;
    pointer[a]++;
  }

  residents1 = (int *)calloc((size_t)npos1 + (size_t)nTotalGrids, sizeof(int));
  heads1 = (int *)calloc((size_t)nTotalGrids, sizeof(int));
  if (residents1 == NULL || heads1 == NULL) {
    goto fail;
  }
  head = 0;
  for (int g = 0; g < nTotalGrids; g++) {
    heads1[g] = head;
    pointer[g] = head;
    head += nResidents1[g] + 1;
    residents1[head - 1] = -1;
  }
  for (int i = 0; i < npos1; i++) {
    int grid[3];
    for (int d = 0; d < 3; d++) {
      grid[d] = frac_to_grid(rpos1[i * 3 + d], ngrid[d]);
    }
    int a = ADDRESS(grid[0], grid[1], grid[2]);
    residents1[pointer[a]] = i;
    pointer[a]++;
  }

  int nGridPairs = gridpairlist(ngrid, False, &gridPairs);
  if (nGridPairs < 0) {
    goto fail;
  }

  int estim = 0;
  for (int i = 0; i < nGridPairs; i++) {
    int g0 = gridPairs[i * 2 + 0];
    int g1 = gridPairs[i * 2 + 1];
    if (accum_prod(&estim, nResidents0[g0], nResidents1[g1]) < 0) {
      goto fail;
    }
  }
  /* NOTE: same-cell pairs are already included in gridPairs when single=False;
     the extra loop below over-allocates (kept for now; cleaned in a later pass). */
  for (int g = 0; g < nTotalGrids; g++) {
    if (accum_prod(&estim, nResidents0[g], nResidents1[g]) < 0) {
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
    for (int j = 0; j < nResidents0[g0]; j++) {
      int r0 = residents0[heads0[g0] + j];
      for (int k = 0; k < nResidents1[g1]; k++) {
        int r1 = residents1[heads1[g1] + k];
        (*pairs)[nPairs * 2 + 0] = r0;
        (*pairs)[nPairs * 2 + 1] = r1;
        nPairs++;
      }
    }
  }

  free(gridPairs);
  free(nResidents0);
  free(nResidents1);
  free(residents0);
  free(residents1);
  free(pointer);
  free(heads0);
  free(heads1);
  return nPairs;

fail:
  free(gridPairs);
  free(nResidents0);
  free(nResidents1);
  free(residents0);
  free(residents1);
  free(pointer);
  free(heads0);
  free(heads1);
  free(*pairs);
  *pairs = NULL;
  return -1;
}


int
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



int
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

