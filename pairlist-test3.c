#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"

//read Gromacs-type data and find the OH covalent bonds.

void
Test(){
  float cell[3];
  float* Hatoms;
  float* Oatoms;
  int nAtoms, nHatoms, nOatoms;
  FILE *file = fopen("pairlist-test3.gro", "r");
  char line[1024];
  fgets( line, sizeof(line), file);
  fgets( line, sizeof(line), file);
  nAtoms = atoi(line);
  nHatoms = 0;
  nOatoms = 0;
  Hatoms = (float*) malloc( sizeof(float) * nAtoms*3 );
  Oatoms = (float*) malloc( sizeof(float) * nAtoms*3 );
  for(int i=0;i<nAtoms;i++){
    fgets(line, sizeof(line), file);
    float x,y,z;
    int   id;
    char  label[14],label2[14];
    sscanf(line, "%s %s %d %f %f %f\n", label, label2, &id, &x, &y, &z);
    //printf("[%s][%d][%f][%f][%f]\n",label2,id,x,y,z);
    if ( strncmp( label2, "HW", 2 ) == 0 ){
      Hatoms[nHatoms*3+0] = x;
      Hatoms[nHatoms*3+1] = y;
      Hatoms[nHatoms*3+2] = z;
      nHatoms ++;
    }
    else if ( strncmp( label2, "OW", 2 ) == 0 ){
      Oatoms[nOatoms*3+0] = x;
      Oatoms[nOatoms*3+1] = y;
      Oatoms[nOatoms*3+2] = z;
      nOatoms ++;
    }
  }
  assert(nOatoms==432);
  assert(nHatoms==864);
  fgets(line, sizeof(line), file);
  sscanf(line, "%f %f %f\n", &cell[0], &cell[1], &cell[2]);
  //printf(line, "%f %f %f\n", cell[0], cell[1], cell[2]);
  int* pairs;
  //find pairs of H and O within 0.12 nm, i.e. covalent OH bonds.
  int nPairs = pairlist(nHatoms, Hatoms, nOatoms, Oatoms, 0.12, cell, &pairs);
  for(int i=0;i<nPairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    float sum2 = 0.0;
    for(int d=0;d<3;d++){
      float delta = Hatoms[r0*3+d] - Oatoms[r1*3+d];
      delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
      sum2 += delta*delta;
    }
    printf("%d %d %f\n", r0,r1, sum2);
  }
  printf("%d nPairs\n", nPairs);
  assert(nPairs==864);
  free(pairs);

  //find pairs of H and O within 0.25 nm, i.e. hydrogen bonds.
  nPairs = pairlist(nHatoms, Hatoms, nOatoms, Oatoms, 0.25, cell, &pairs);
  int nHBs = 0;
  for(int i=0;i<nPairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    float sum2 = 0.0;
    for(int d=0;d<3;d++){
      float delta = Hatoms[r0*3+d] - Oatoms[r1*3+d];
      delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
      sum2 += delta*delta;
    }
    //if the pair is not covalent,
    if ( 0.12*0.12 < sum2 ){
      printf("%d %d %f\n", r0,r1, sum2);
      nHBs ++;
    }
  }
  printf("%d nHBs\n", nHBs);
  assert(nHBs==864);
  free(pairs);



  free(Hatoms);
  free(Oatoms);
}



int main(int argc, char * argv[])
{
  Test();
}


