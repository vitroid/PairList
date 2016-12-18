#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define ADDRESS(x,y,z) (((z)*GY + (y))*GX + (x))

int
gridpairlist(int ngrid[3], int **pairs)
/* 
given:
    ngrid: number of grids in three dimensions
returns:
    return value: number of grid pairs
    pairs: newly allocated list of grid pairs
*/
{
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  const int nTotalGrids = GX*GY*GZ;
  int npair = 0;
  (*pairs) = (int*) malloc(sizeof(int)*14*nTotalGrids*2);
  //make neighboring grid list
  for(int ix=0;ix<GX;ix++){
    for(int iy=0;iy<GY;iy++){
      for(int iz=0;iz<GZ;iz++){
        int a = ADDRESS(ix,iy,iz); //center
        for(int kx=ix-1;kx<=ix+1;kx++){
          int jx = (kx+GX) % GX;
          for(int ky=iy-1;ky<=iy+1;ky++){
            int jy = (ky+GY) % GY;
            for(int kz=iz-1;kz<=iz+1;kz++){
              int jz = (kz+GZ) % GZ;
              int b = ADDRESS(jx,jy,jz);
              if ( a<b ){
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


int
pairlist(int nAtoms, float *atoms, float threshold, float cell[3], int **pairs)
/* 
given:
    nAtoms: number of atoms
    atoms: positions of atoms
    threshold: two atoms are considered to be paired when their distance is below it.
    cell: rectangular cell shape
returns:
    return value: number of pairs
    pairs: newly allocated list of pairs
*/
{
  int test = 0;
  //determine the grid size
  int ngrid[3];
  for(int d=0;d<3;d++){
    ngrid[d] = (int) floor(cell[d] / threshold);
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  const int nTotalGrids = GX*GY*GZ;
  
  //count the number of residents in each grid cells.
  int nResidents[nTotalGrids];
  for(int i=0;i<nTotalGrids;i++) nResidents[i] = 0;
  for(int i=0;i<nAtoms;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms[i*3+d] / cell[d];
      r[d] -= floor(r[d]); //Periodic boundary condition
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    nResidents[ADDRESS(grid[0],grid[1],grid[2])] ++;
  }

  //make the resident list
  //resident list is serialized to reduce the memory usage
  //head pointer to the list of residents in a grid cell is in headOfList.
  //atom ids are put in the residents[] and is terminated with -1.
  int residents[nAtoms + nTotalGrids];
  int headOfList[nTotalGrids];
  int pointer[nTotalGrids];
  int head = 0;
  for(int g=0;g<nTotalGrids;g++){
    headOfList[g] = head;
    pointer[g] = head;
    head += nResidents[g]+1;
    residents[head-1] = -1;
  }
  for(int i=0;i<nAtoms;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms[i*3+d] / cell[d];
      r[d] -= floor(r[d]);
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    int a = ADDRESS(grid[0],grid[1],grid[2]);
    residents[pointer[a]] = i;
    pointer[a] ++;
  }
  if(test){
    //just print
    for(int g=0;g<nTotalGrids;g++){
      printf("Grid %d\n", g);
      printf("Number of residents: %d\n", nResidents[g]);
      printf("First resident: %d\n", residents[headOfList[g]]);
      printf("Last resident: %d\n", residents[headOfList[g] + nResidents[g] - 1]);
      printf("Terminator (must be -1): %d\n", residents[headOfList[g] + nResidents[g]]);
    }
  }
  //make the neighboring grid pair list
  int* gridPairs;
  int nGridPairs = gridpairlist(ngrid, &gridPairs);
  //rough estimate of the number of atom pairs
  int rough = 0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    int g1 = gridPairs[i*2+1];
    rough += nResidents[g0] * nResidents[g1];
  }
  for(int g=0;g<nTotalGrids;g++){
    rough += nResidents[g]*nResidents[g] / 2;
  }

  *pairs = (int*) malloc(sizeof(int) * rough*2);
  int nPairs=0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    for(int j=0; j<nResidents[g0]; j++){
      int r0 = residents[headOfList[g0] + j];
      int g1 = gridPairs[i*2+1];
      for(int k=0; k<nResidents[g1]; k++){
        int r1 = residents[headOfList[g1] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms[r0*3+d] - atoms[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if (sum2 < threshold*threshold){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  for(int g=0;g<nTotalGrids;g++){
    for(int j=0; j<nResidents[g]; j++){
      int r0 = residents[headOfList[g] + j];
      for(int k=j+1; k<nResidents[g]; k++){
        int r1 = residents[headOfList[g] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms[r0*3+d] - atoms[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if (sum2 < threshold*threshold){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  //deallocate allocated memories
  free(gridPairs);
  return nPairs;
}
  
