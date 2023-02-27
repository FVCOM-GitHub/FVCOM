!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

SUBROUTINE LOAD_GRID
!===============================================================================!
! GET THE GRID CONNECTIVITY AND OBC NODE LIST FOR DOMAIN
!    DECOMPOSITION AND GENMAP
!===============================================================================!
  USE CONTROL
  USE MOD_INPUT
  USE MOD_OBCS
  USE ALL_VARS
  USE MOD_NESTING
  IMPLICIT NONE
  integer i


  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE("hotstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR HOTSTART                        !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL LOAD_RESTART_GRID(NVG) ! SET DIMENSIONS HERE
          
     CALL LOAD_RESTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
        
     CALL LOAD_RESTART_LSF_GRID(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)


     !=================================================
     ! CRASHSTART 
  CASE("crashrestart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR CRASHSTART                      !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL LOAD_RESTART_GRID(NVG)  ! SET DIMENSIONS HERE
 
     CALL LOAD_RESTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
     
     CALL LOAD_RESTART_LSF_GRID(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)  

     !=================================================
     ! COLDSTART 
  CASE("coldstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     ! SET DIMENSIONS HERE
     CALL LOAD_COLDSTART_SIGMA 

     CALL LOAD_COLDSTART_GRID(NVG)
     
     CALL LOAD_COLDSTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)

     CALL LOAD_COLDSTART_LSF(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)

  END SELECT
  

END SUBROUTINE LOAD_GRID
