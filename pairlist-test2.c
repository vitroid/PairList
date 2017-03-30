#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"

//read NX4A-type data and find the HB pairs.

void
Test2(){
  float cell[3];
  float* atoms;
  int nAtoms;
  FILE *file = fopen("pairlist-test2.data", "r");
  char line[1024];
  while( fgets( line, sizeof(line), file) ){
    if (strncmp(line, "@BOX3", 5) == 0 ){
      fgets(line, sizeof(line), file);
      sscanf(line, "%f %f %f\n", &cell[0], &cell[1], &cell[2]);
    }
    else if ( strncmp(line, "@NX4A", 5 ) == 0){
      fgets(line, sizeof(line), file);
      nAtoms = atoi(line);
      atoms = (float*) malloc( sizeof(float) * nAtoms*3 );
      for(int i=0;i<nAtoms;i++){
        fgets(line, sizeof(line), file);
        sscanf(line, "%f %f %f\n", &atoms[i*3+0], &atoms[i*3+1], &atoms[i*3+2]);
      }
      break;
    }
  }
  int* pairs;
  int nPairs = pairlist(nAtoms, atoms, 0, NULL, 0.0, 3.0, cell, &pairs);
  for(int i=0;i<nPairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    float sum2 = 0.0;
    for(int d=0;d<3;d++){
      float delta = atoms[r0*3+d] - atoms[r1*3+d];
      delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
      sum2 += delta*delta;
    }
    printf("%d %d %f\n", r0,r1, sum2);
  }
  assert(nPairs==864);
  free(pairs);
  free(atoms);
}



int main(int argc, char * argv[])
{
  Test2();
}


