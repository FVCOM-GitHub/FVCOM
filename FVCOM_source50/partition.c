
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * 
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partdmesh.c,v 1.1 1998/11/27 17:59:38 karypis Exp $
 *
 */

/*
 * in frename.c, 6 subroutine names are mapped to the following suroutine
 * this supports the situation where the fortran compiler namespace is different
 * which results in the fortran routine calling a variety of different names
 * for Partition including:
 * 	partition
 * 	PARTITION
 * 	partition_
 * 	PARTITION_
 * 	partition__
 * 	PARTITION__
 * if the fortran compiler gives an error during linking that the routine
 * partition can't be found and you are sure you are linking in the metis 
 * libraries correctly, check carefully the name of the subroutine the compiler
 * is looking for.  If it is a different variant on "partition", add this 
 * variant at the bottom of "frename.c" and recompile the library*/

/* Adapted from partdmesh.c by G. Cowles (April, 2016) for the purpose of
 * binding FVCOM (fortran) and METIS (c) 
 */ 

#include <stdlib.h>
#include <metis.h>
#include <stdio.h>

#define LTERM                   (void **) 0

void partition_(int* npartsin, int* nein, int* nnin, int* ncommonin,  int* nv[], int* optionsin[], int* epartout[], int* objvalout)
{
  int i, j, nparts, edgecut, nn, ne, ncommon, objval;
  idx_t *eptr, *eind, *epart, *npart, *options;
   
  ncommon = *ncommonin;
  nparts = *npartsin; 
  nn = *nnin;
  ne = *nein;

  /* local vars */
  nn = *nnin;
  ne = *nein;
  eind = (malloc(sizeof(int)*(3*ne)));
  eptr = (malloc(sizeof(int)*(ne+1)));
  epart = (malloc(sizeof(int)*(ne)));
  npart = (malloc(sizeof(int)*(nn)));
  options = (malloc(sizeof(int)*(40)));

  for (j=3*ne, i=0; i<j; i++) {
      eind[i] = (*nv)[i];
  }

  for (j=ne+1, i=0; i<j; i++){
      eptr[i] = i*3+1; 
  }

  for (j=40, i=0; i < j; i++) {
      options[i] = (*optionsin)[i];
  }


  METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, options, &objval, epart, npart);
  printf("edgecut %d\n",objval); 

  // set return values
  (*objvalout) = objval;

  for (j=ne, i=0; i<j; i++) {
       (*epartout)[i] = epart[i];
  }

  // free up the memory
  free(epart);
  free(eptr);
  free(eind);
  free(npart);
  free(options);

}

