/*==========================================================================
* Copyright (c) 2007, The University of Massachusetts Dartmouth 
* Produced at the School of Marine Science & Technology 
* Marine Ecosystem Dynamics Modeling group
* All rights reserved.
*
* FVCOM has been developed by the joint UMASSD-WHOI research team. For 
* details of authorship and attribution of credit please see the FVCOM
* technical manual or contact the MEDM group.
*
* 
* This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
* The full copyright notice is contained in the file COPYRIGHT located in the 
* root directory of the FVCOM code. This original header must be maintained
* in all distributed versions.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  
*
*---------------------------------------------------------------------------
* CVS VERSION INFORMATION
* $Id$
* $Name$
* $Revision$
*===========================================================================*/

#include <stdio.h>

#define maxfuncs 20

/* Can you pass fortran pointers through c? */

typedef struct {
  void (*func)();
} VOIDFUNC_TABLE;


VOIDFUNC_TABLE voidfuncs[maxfuncs];

static int nvoidfuncs =0;

void register_func(void (*func)(), int *i, int *stat) {
   *stat = -1;
     if(nvoidfuncs <= maxfuncs)
       {
         /* Convert from fortran indexing */
         voidfuncs[nvoidfuncs++].func = func; /*Insert at nvoidfuncs */
         *i = nvoidfuncs; /*Return value nvoidfuncs+1 */
         *stat =0;
       }
}

void register_func_(void (*func)(), int *i, int *stat) {
   *stat = -1;
     if(nvoidfuncs <= maxfuncs)
       {
         /* Convert from fortran indexing */
         voidfuncs[nvoidfuncs++].func = func; /*Insert at nvoidfuncs */
         *i = nvoidfuncs; /*Return value nvoidfuncs+1 */
         *stat =0;
       }
}



void call_func(int *i, int *stat){
  *stat = -1;
  int j = (*i)-1; /* Convert from fortran indexing */
  if ( j >= 0 && j < nvoidfuncs)
    {
      voidfuncs[j].func();
      *stat=0;
    }
}


void call_func_(int *i, int *stat){
  *stat = -1;
  int j = (*i)-1; /* Convert from fortran indexing */
  if ( j >= 0 && j < nvoidfuncs)
    {
      voidfuncs[j].func();
      *stat=0;
    }
}
