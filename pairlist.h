#ifndef __PAIRLIST_H__
#define __PAIRLIST_H__
int pairlist(int nAtoms, float *atoms, float threshold, float cell[3], int **pairs);
int pairlist2(int nAtoms0, float *atoms0, int nAtoms1, float *atoms1, float threshold, float cell[3], int **pairs);


#endif
