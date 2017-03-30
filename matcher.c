#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"

//read Gromacs-type data and find the OH covalent bonds.

int
LoadGRO(float** Oatoms, float* cell){
  int nAtoms, nOatoms;
  FILE *file = fopen("TPPI_300K.gro", "r");
  char line[1024];
  fgets( line, sizeof(line), file);
  fgets( line, sizeof(line), file);
  nAtoms = atoi(line);
  nOatoms = 0;
  *Oatoms = (float*) malloc( sizeof(float) * nAtoms*3 );
  for(int i=0;i<nAtoms;i++){
    fgets(line, sizeof(line), file);
    float x,y,z;
    int   id;
    char  label[99],label2[99];
    sscanf(line, "%10s%5s%5d%8f%8f%8f\n", label, label2, &id, &x, &y, &z);
    if ( strncmp( label2, "Ow", 2 ) == 0 ){
      (*Oatoms)[nOatoms*3+0] = x;
      (*Oatoms)[nOatoms*3+1] = y;
      (*Oatoms)[nOatoms*3+2] = z;
      nOatoms ++;
    }
  }
  //assert(nOatoms==nAtoms/3);
  fgets(line, sizeof(line), file);
  sscanf(line, "%f %f %f\n", &cell[0], &cell[1], &cell[2]);
  return nOatoms;
}


void
MakeNeighborList(int natoms, int npairs, int* pairs, 
		 //return values
		 int* nnei, int** nei)
{
  //make neighborlist
  for(int i=0;i<natoms; i++){
    nnei[i] = 0;
  }
  for(int i=0;i<npairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    nnei[r0] += 1;
    nnei[r1] += 1;
  }
  for(int i=0;i<natoms; i++){
    nei[i] = (int*)malloc(nnei[i]*sizeof(int));
    nnei[i] = 0;
  }
  for(int i=0;i<npairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    nei[r0][nnei[r0]] = r1;
    nnei[r0]+=1;
    nei[r1][nnei[r1]] = r0;
    nnei[r1]+=1;
  }
}


void Test()
{
  float cell[3];
  float *Oatoms;
  int nOatoms = LoadGRO(&Oatoms, cell);
  printf("%d %f %f %f\n", nOatoms, cell[0], cell[1], cell[2]);
  //Find the triangle of a and b cell vectors
  float a = 0.939*sqrt(2);
  //float b = 0.939*sqrt(2);  same as a
  float c = 0.738;
  float ab = 0.939*2;
  float ac = sqrt(a*a + c*c);
  int* pairsA;
  int* pairsC;
  int* pairsAB;
  int* pairsAC;
  int nPairsA = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, a*0.97, a*1.03, cell, &pairsA);
  int nPairsC = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, c*0.97, c*1.03, cell, &pairsC);
  int nPairsAB = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, ab*0.97, ab*1.03, cell, &pairsAB);
  int nPairsAC = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, ac*0.97, ac*1.03, cell, &pairsAC);
  int nneiA[nOatoms];
  int* neiA[nOatoms];
  int nneiC[nOatoms];
  int* neiC[nOatoms];
  int nneiAB[nOatoms];
  int* neiAB[nOatoms];
  int nneiAC[nOatoms];
  int* neiAC[nOatoms];
  MakeNeighborList(nOatoms, nPairsA, pairsA, nneiA, neiA);
  MakeNeighborList(nOatoms, nPairsC, pairsC, nneiC, neiC);
  MakeNeighborList(nOatoms, nPairsAB, pairsAB, nneiAB, neiAB);
  MakeNeighborList(nOatoms, nPairsAC, pairsAC, nneiAC, neiAC);
  //find triangle PQR that matches the shape
  int ntri=0;
  for(int p=0; p<nOatoms; p++){
    for(int i=0; i<nneiA[p]; i++){
      int q = neiA[p][i];
      for(int j=0; j<nneiA[p]; j++){
	int r = neiA[p][j];
	int found = 0;
	for(int k=0; k<nneiAB[q]; k++){
	  if ( neiAB[q][k] == r ){
	    found = 1;
	    break;
	  }
	}
	if ( found ){
	  //printf("%d %d %d\n", p,q,r);
	  ntri += 1;
	}
      }
    }
  }
  printf("%d ntri\n", ntri);
  /*
  for(int i=0;i<nPairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    float sum2 = 0.0;
    for(int d=0;d<3;d++){
      float delta = Oatoms[r0*3+d] - Oatoms[r1*3+d];
      delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
      sum2 += delta*delta;
    }
    printf("%d %d %f\n", r0,r1, sqrt(sum2));
  }
  */
  printf("%d nPairs A\n", nPairsA);
  printf("%d nPairs C\n", nPairsC);
  printf("%d nPairs AB\n", nPairsAB);
  printf("%d nPairs AC\n", nPairsAC);
  //assert(nPairs==864);
  free(pairsA);
  free(pairsC);
  free(pairsAB);
  free(pairsAC);

  free(Oatoms);
}



int main(int argc, char * argv[])
{
  Test();
}


