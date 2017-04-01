#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"
#include "bst.h"

//read Gromacs-type data and find the OH covalent bonds.

//benchmark on vitroid-black 2017-3-31
//original matcher (no bst):             user	1m28.847s
//matcher with bst:                      user	1m48.723s
//matcher with bst (size() is improved): user	1m48.942s

#define max(A,B) ((A)>(B)?(A):(B))
int
isClose(float x, float y)
{
  float e = (x-y)/max(x,y);
  return (-1e-3 < e) && (e < +1e-3);
}




void
yaplot_linerel(float* a, float* d)
{
  printf("l %f %f %f %f %f %f\n", a[0],a[1],a[2], a[0]+d[0],a[1]+d[1], a[2]+d[2]);
}

int
LoadGRO(FILE* file, float** Oatoms, float* cell){
  int nAtoms, nOatoms;
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


float dot(float* x, float* y)
{
  float sum =0;
  for(int d=0;d<3;d++){
    sum += x[d]*y[d];
  }
  return sum;
}

  
float vector_length(float* v)
{
  return sqrt(dot(v,v));
}

void
cross(float*x, float* y, float* z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}


void
sub(float*x, float* y, float* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] -y[d];
  }
}


void
mul(float*x, float a, float* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] * a;
  }
}


void
regularize(float* x, float* y, float* z)
{
  float z0[3];
  for(int d=0;d<3;d++)
    z0[d] = z[d];
  float error = dot(x,y);
  float ye[3];
  mul(y,error/2,ye);
  float x_ort[3];
  sub(x,ye,x_ort);
  float xe[3];
  mul(x,error/2,xe);
  float y_ort[3];
  sub(y,xe,y_ort);
  float z_ort[3];
  cross(x_ort, y_ort, z_ort);
  mul(x_ort, 0.5*(3.0-dot(x_ort,x_ort)), x);
  mul(y_ort, 0.5*(3.0-dot(y_ort,y_ort)), y);
  mul(z_ort, 0.5*(3.0-dot(z_ort,z_ort)), z);
  //allow mirroring
  if (dot(z0, z) < 0){
    mul(z, -1, z);
  }
}


float
det(float* A, float* B, float* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  float a = A[0];
  float b = B[0];
  float c = C[0];
  float d = A[1];
  float e = B[1];
  float f = C[1];
  float g = A[2];
  float h = B[2];
  float k = C[2];
  return a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
}





void
inv(float* A, float* B, float* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  float D = det(A,B,C);
  float a = A[0];
  float b = B[0];
  float c = C[0];
  float d = A[1];
  float e = B[1];
  float f = C[1];
  float g = A[2];
  float h = B[2];
  float k = C[2];
  A[0] =  (e * k - f * h) / D;
  B[0] = -(b * k - c * h) / D;
  C[0] =  (b * f - c * e) / D;
  A[1] = -(d * k - f * g) / D;
  B[1] =  (a * k - c * g) / D;
  C[1] = -(a * f - c * d) / D;
  A[2] =  (d * h - e * g) / D;
  B[2] = -(a * h - b * g) / D;
  C[2] =  (a * e - b * d) / D;
}

    





int
LoadAR3R(FILE* file, float** atoms, float* cell){
  int nAtoms;
  char line[1024];
  while ( NULL != fgets( line, sizeof(line), file) ){
    if ( 0 == strncmp(line, "@BOX3", 5) ){
      fgets( line, sizeof(line), file);
      sscanf(line, "%f %f %f\n", &cell[0], &cell[1], &cell[2]);
      cell[0] /= 10; //in nm, for gromacs
      cell[1] /= 10;
      cell[2] /= 10;
    }
    else if ( 0 == strncmp(line, "@AR3R", 5) ){
      fgets( line, sizeof(line), file);
      nAtoms = atoi(line);
      *atoms = (float*) malloc( sizeof(float) * nAtoms*3 );
      float dmin = 1e99;
      int   imin = -1;
      for(int i=0;i<nAtoms;i++){
	fgets(line, sizeof(line), file);
	float x,y,z;
	sscanf(line, "%f %f %f", &x, &y, &z);
	//atoms distribute around the origin.
	x -= floor(x+0.5);
	y -= floor(y+0.5);
	z -= floor(z+0.5);
	float rr = x*x + y*y + z*z;
	if ( rr < dmin ){
	  dmin = rr;
	  imin = i;
	}
	(*atoms)[i*3+0] = x;
	(*atoms)[i*3+1] = y;
	(*atoms)[i*3+2] = z;
      }
      //there must be a molecule at the origin
      for(int i=0;i<nAtoms;i++){
	(*atoms)[i*3+0] -= (*atoms)[imin*3+0];
	(*atoms)[i*3+1] -= (*atoms)[imin*3+1];
	(*atoms)[i*3+2] -= (*atoms)[imin*3+2];
      }
    }
  }
  return nAtoms;
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





int main(int argc, char* argv[])
{
  //usage: matcher grofile error msdmax
  //typical values: error = 0.02, msdmax == 0.9 for TPPI
  if ( argc != 5 ){
    fprintf(stderr, "usage: %s grofile error msdmax\n", argv[0]);
    exit(1);
  }
  float cell[3];
  float *Oatoms;
  FILE *file = fopen(argv[1], "r");
  float err, msdmax;
  sscanf(argv[3], "%f", &err);
  sscanf(argv[4], "%f", &msdmax);
  int nOatoms = LoadGRO(file, &Oatoms, cell);
  fclose(file);

  float unitcell[3];
  float unitcelli[3];
  float *unitatoms; //relative
  file = fopen(argv[2], "r");
  int nunitatoms = LoadAR3R(file, &unitatoms, unitcell);
  fclose(file);
  assert(isClose(unitcell[0], unitcell[1]));
  float radius = vector_length(unitcell)/2;
  for(int d=0;d<3;d++){
    unitcelli[d] = 1.0/unitcell[d];
  }
  fprintf(stderr,"System   %d %f %f %f\n", nOatoms, cell[0], cell[1], cell[2]);
  fprintf(stderr,"Template %d %f %f %f\n", nunitatoms, unitcell[0], unitcell[1], unitcell[2]);
  //Find the tetrahedra of a, b, and c cell vectors
  float a = unitcell[0];
  //float b   same as a
  float c = unitcell[2];
  float ab = sqrt(a*a*2);
  float ac = sqrt(a*a + c*c);
  int* prox;
  int* pairsA;
  int* pairsC;
  int* pairsAB;
  int* pairsAC;
  //atoms of the proximity
  fprintf(stderr, "Preparing the neighbor lists.");
  int nProx   = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, 0.0, radius*(1+err), cell, &prox);
  //pair lists
  fprintf(stderr, ".");
  int nPairsA = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, a*(1-err), a*(1+err), cell, &pairsA);
  fprintf(stderr, ".");
  int nPairsC = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, c*(1-err), c*(1+err), cell, &pairsC);
  fprintf(stderr, ".");
  int nPairsAB = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, ab*(1-err), ab*(1+err), cell, &pairsAB);
  fprintf(stderr, ".");
  int nPairsAC = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, ac*(1-err), ac*(1+err), cell, &pairsAC);
  bnode* nei[nOatoms];
  bnode* neiA[nOatoms];
  bnode* neiC[nOatoms];
  bnode* neiAB[nOatoms];
  bnode* neiAC[nOatoms];
  fprintf(stderr, ":");
  MakeNeighborList(nOatoms, nProx, prox, nei);
  fprintf(stderr, ":");
  MakeNeighborList(nOatoms, nPairsA, pairsA, neiA);
  fprintf(stderr, ":");
  MakeNeighborList(nOatoms, nPairsC, pairsC, neiC);
  fprintf(stderr, ":");
  MakeNeighborList(nOatoms, nPairsAB, pairsAB, neiAB);
  fprintf(stderr, ":");
  MakeNeighborList(nOatoms, nPairsAC, pairsAC, neiAC);
  free(prox);
  free(pairsA);
  free(pairsC);
  free(pairsAB);
  free(pairsAC);
  fprintf(stderr, "Done.\n");
  fprintf(stderr,"%d nProx A\n", nProx);
  fprintf(stderr,"%d nPairs A\n", nPairsA);
  fprintf(stderr,"%d nPairs C\n", nPairsC);
  fprintf(stderr,"%d nPairs AB\n", nPairsAB);
  fprintf(stderr,"%d nPairs AC\n", nPairsAC);
  //find triangle PQR that matches the shape
  int ntet=0;
  for(int p=0; p<nOatoms; p++){
    fprintf(stderr, "\r%.1f %%", p*100./nOatoms);
    int nnA = size(neiA[p]);
    //fprintf(stderr,"size of neiA is %d\n", nnA);
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
	      //fprintf(stderr,"%d %d %d\n", p,q,r);

	      float dq[3], dr[3], ds[3];
	      for(int d=0;d<3;d++){
		dq[d] = Oatoms[q*3+d] - Oatoms[p*3+d];
		dq[d] -= floor( dq[d] / cell[d] + 0.5 ) * cell[d];
		dr[d] = Oatoms[r*3+d] - Oatoms[p*3+d];
		dr[d] -= floor( dr[d] / cell[d] + 0.5 ) * cell[d];
		ds[d] = Oatoms[s*3+d] - Oatoms[p*3+d];
		ds[d] -= floor( ds[d] / cell[d] + 0.5 ) * cell[d];
	      }
	      /*
	      yaplot_linerel(&Oatoms[p*3], dq);
	      yaplot_linerel(&Oatoms[p*3], dr);
	      yaplot_linerel(&Oatoms[p*3], ds);
	      */
	      //mathcing without alignment.
	      //find the rotation matrix from 
	      //A is the unit cell
	      //B is the unit cell in the atomic configuration
	      //R is the rotation for A
	      //Assume A and B has common origin. i.e. no translation
	      //Then, if both A and B are not distorted,
	      //B = RA, i.e. R = B A^-1
	      //if A and B are stacks of row vectors,
	      //B = AR, i.e. R = A^-1 B
	      //In reality, R wont be a regular rotation matrix, so we have to regularize it.
	      //Here is a simple idea.
	      //http://stackoverflow.com/questions/23080791/eigen-re-orthogonalization-of-rotation-matrix
	      //In the present case, unit cell is orthogonal, so inverse is trivial.
	      //Assume A and B are row vectors
	      float Rx[3], Ry[3], Rz[3];
	      mul(dq, unitcelli[0], Rx);
	      mul(dr, unitcelli[1], Ry);
	      mul(ds, unitcelli[2], Rz);
	      //printf("x0 %f %f %f %f\n", Rx[0],Rx[1],Rx[2], dot(Rx,Rx));
	      //printf("y0 %f %f %f %f\n", Ry[0],Ry[1],Ry[2], dot(Ry,Ry));
	      //printf("z0 %f %f %f %f\n%f\n", Rz[0],Rz[1],Rz[2], dot(Rz,Rz), det(Rx,Ry,Rz));
	      regularize(Rx,Ry,Rz);
	      //printf("x1 %f %f %f %f\n", Rx[0],Rx[1],Rx[2], dot(Rx,Rx));
	      //printf("y1 %f %f %f %f\n", Ry[0],Ry[1],Ry[2], dot(Ry,Ry));
	      //printf("z1 %f %f %f %f\n%f\n\n", Rz[0],Rz[1],Rz[2],dot(Rz,Rz), det(Rx,Ry,Rz));
	      //We get R' by regularization, and B' = A R'
	      float ARx[3], ARy[3], ARz[3];
	      mul(Rx, unitcell[0], ARx);
	      mul(Ry, unitcell[1], ARy);
	      mul(Rz, unitcell[2], ARz);
	      //this does not include p itself!
	      int* neighbors = get_array(nei[p]);
	      neighbors[size(nei[p])] = p; // add itself in place of the useless terminator
	      float msd = 0.0;
	      int partners[nunitatoms];
	      for(int l=0; l<nunitatoms; l++){
		float u[3];
		u[0] = dot(&unitatoms[l*3], ARx) + Oatoms[p*3+0];
		u[1] = dot(&unitatoms[l*3], ARy) + Oatoms[p*3+1];
		u[2] = dot(&unitatoms[l*3], ARz) + Oatoms[p*3+2];
		float dmin = 1e99;
		for(int m=0; m<size(nei[p])+1; m++){
		  float dd[3];
		  int ne = neighbors[m];
		  sub(u, &Oatoms[ne*3], dd);
		  float L = dot(dd,dd);
		  if ( L < dmin ){
		    dmin = L;
		    partners[l] = ne;
		  }
		}
		msd += dmin;
	      }
	      if ( msd < msdmax ){
		/*
		for(int l=0; l<nunitatoms; l++){
		  float u[3];
		  u[0] = dot(&unitatoms[l*3], ARx) + Oatoms[p*3+0];
		  u[1] = dot(&unitatoms[l*3], ARy) + Oatoms[p*3+1];
		  u[2] = dot(&unitatoms[l*3], ARz) + Oatoms[p*3+2];
		  printf("r 0.01\n");
		  printf("@ 3\n");
		  printf("c %f %f %f\n", u[0],u[1],u[2]);
		  int pa = partners[l];
		  printf("c %f %f %f\n", Oatoms[pa*3+0],Oatoms[pa*3+1],Oatoms[pa*3+2]);
		  printf("l %f %f %f %f %f %f\n", u[0],u[1],u[2], Oatoms[pa*3+0],Oatoms[pa*3+1],Oatoms[pa*3+2]);
		}
		exit(0);
		*/
		printf("%d ", nunitatoms);
		for(int l=0;l<nunitatoms;l++){
		  printf("%d ", partners[l]);
		}
		printf("%f %d ", msd, p);//error and the origin atom
		printf("%f %f %f ", Rx[0],Rx[1],Rx[2]);
		printf("%f %f %f ", Ry[0],Ry[1],Ry[2]);
		printf("%f %f %f\n", Rz[0],Rz[1],Rz[2]);

		
	      }
	      free(neighbors);
		  
	      ntet += 1;
	    }
	  }
	  free(nC);
	}
      }
    }
    free(nA);
  }
  //dispose memory as soon as possible
  for(int i=0;i<nOatoms; i++){
    dispose(nei[i]);
    dispose(neiA[i]);
    dispose(neiC[i]);
    dispose(neiAB[i]);
    dispose(neiAC[i]);
  }
  fprintf(stderr,"%d ntet\n", ntet);
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
    fprintf(stderr,"%d %d %f\n", r0,r1, sqrt(sum2));
  }
  */
  //assert(nPairs==864);

  free(Oatoms);
}




