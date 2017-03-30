#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"
#include "bst.h"

//read Gromacs-type data and find the OH covalent bonds.

//benchmark on vitroid-black 2017-3-31
//original matcher (no bst): user	1m28.847s
//matcher with bst: user	1m48.723s

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
		 bnode** nei)
{
  //make neighborlist
  for(int i=0;i<natoms; i++){
    nei[i] = NULL;
  }
  for(int i=0;i<npairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    nei[r0] = insert(nei[r0], r1);
    nei[r1] = insert(nei[r1], r0);
  }
}


int
isin(int value, int size, int* arr)
{
  for(int i=0;i<size; i++){
    if (arr[i] == value){
      return 1;
    }
  }
  return 0;
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
  bnode* neiA[nOatoms];
  bnode* neiC[nOatoms];
  bnode* neiAB[nOatoms];
  bnode* neiAC[nOatoms];
  MakeNeighborList(nOatoms, nPairsA, pairsA, neiA);
  MakeNeighborList(nOatoms, nPairsC, pairsC, neiC);
  MakeNeighborList(nOatoms, nPairsAB, pairsAB, neiAB);
  MakeNeighborList(nOatoms, nPairsAC, pairsAC, neiAC);
  //find triangle PQR that matches the shape
  int ntet=0;
  for(int p=0; p<nOatoms; p++){
    int nnA = size(neiA[p]);
    //printf("size of neiA is %d\n", nnA);
    int* nA = get_array(neiA[p]);
    for(int i=0; i<nnA; i++){
      int q = nA[i];
      for(int j=0; j<nnA; j++){
	int r = nA[j];
	if ( lookup(neiAB[q], r) ){
	  
	  int* nC = get_array(neiC[p]);
	  int nnC = size(neiC[p]);
	  for(int k=0; k<nnC; k++){
	    int s = nC[k];
	    if ( lookup(neiAC[q], s) && lookup(neiAC[r], s) ){
	      //printf("%d %d %d\n", p,q,r);
	      ntet += 1;
	    }
	  }
	  free(nC);
	}
      }
    }
    free(nA);
  }
  for(int i=0;i<nOatoms; i++){
    dispose(neiA[i]);
    dispose(neiC[i]);
    dispose(neiAB[i]);
    dispose(neiAC[i]);
  }
  printf("%d ntet\n", ntet);
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


