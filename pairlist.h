#ifndef __PAIRLIST_H__
#define __PAIRLIST_H__
int pairlist(int nAtoms0, float *atoms0, int nAtoms1, float *atoms1, float lower, float higher, float cell[3], int **pairs);
int Pairs(int npos, float* rpos, int ngrid[3], int **pairs);
int Pairs2(int npos0, float *rpos0, int npos1, float *rpos1, int ngrid[3], int **pairs);


#endif
