#ifndef __PAIRLIST_H__
#define __PAIRLIST_H__
#define PL_FLOAT double
int pairlist(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);
/* Pairs/Pairs2: ngrid[d] must be >= 1.
   On success return pair count (>=0); *pairs is malloc'd (or NULL if count==0).
   On error return -1 and set *pairs to NULL. */
int Pairs(int npos, PL_FLOAT* rpos, int ngrid[3], int **pairs);
int Pairs2(int npos0, PL_FLOAT *rpos0, int npos1, PL_FLOAT *rpos1, int ngrid[3], int **pairs);


#endif
