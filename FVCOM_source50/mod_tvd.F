! -----------------------------------------------------------------------------
! The TVD option for FVCOM was developed by Akvaplan-niva in 2018.
! Preliminary results suggest that TVD advection achieves less numerical
! diffusion and better shape preserving properties than the standard scheme. 
!
! This scheme is also 2nd order accurate, the difference to the standard
! scheme lies in how we interpolate values to the control-volume walls.
!
! We recommend those interested in a deeper understanding of TVD schemes can be
! applied for unstructured mesh models to have a look at Zhang et. al. (2015), 
! Zhang et. al. (2016) and Darwish and Moukalled (2003) 
!
! Zhang, Di, et al. "A refined r‐factor algorithm for TVD schemes on arbitrary 
! unstructured meshes." 
! International Journal for Numerical Methods in Fluids 80.2 (2016): 105-139.
!
! Zhang, Zhuo, et al. "A new r‐ratio formulation for TVD schemes for 
! vertex‐centered FVM on an unstructured mesh." 
! International Journal for Numerical Methods in Fluids 81.12 (2016): 741-764.
!
! Darwish, M. S., and F. Moukalled. "TVD schemes for unstructured grids." 
! International Journal of heat and mass transfer 46.4 (2003): 599-611.
!
! - Håvard Espenes, Ole Anders Nøst
! -----------------------------------------------------------------------------

MODULE MOD_TVD

   USE MOD_PREC
   IMPLICIT NONE
   SAVE
!
!--Parameters for TVD advection                 
!
   REAL(SP), ALLOCATABLE :: DELF(:)  

   REAL(SP), ALLOCATABLE :: Anear_node(:), Bnear_node(:)
   REAL(SP), ALLOCATABLE :: YUAdist(:), YUBdist(:), XUAdist(:), XUBdist(:)

CONTAINS

!==============================================================================!
  SUBROUTINE SETUP_TVD
    USE ALL_VARS
    USE LIMS
    IMPLICIT NONE
    REAL(SP), DIMENSION(1:M)             :: dstA,dstB
    REAL(SP), DIMENSION(1:NCV)        :: XUA,XUB,YUA,YUB
    REAL(SP) :: LX,LY,minA,minB
    INTEGER  :: I,I2,IA,IB

    ALLOCATE(Anear_node(1:NCV))   ; Anear_node     = ZERO
    ALLOCATE(Bnear_node(1:NCV))   ; Bnear_node     = ZERO
    ALLOCATE(XUAdist(1:NCV))      ; XUAdist    = ZERO
    ALLOCATE(YUAdist(1:NCV))      ; YUAdist    = ZERO
    ALLOCATE(XUBdist(1:NCV))      ; XUBdist    = ZERO
    ALLOCATE(YUBdist(1:NCV))      ; YUBdist    = ZERO



  ! --------------------- Calculate TVD parameters ----------------------------- !
  DO I=1,NCV
      IA = NIEC(I,1) ! The node to the left of the controllvolume-edge
      IB = NIEC(I,2) ! The node to the right of the controllvolume-edge

      LX = VX(IB)-VX(IA) ! X-distance between IA and IB
      LY = VY(IB)-VY(IA) ! Y-distance between IA and IB

  ! If A is upstream
      XUA(I) = VX(IA)-LX ! x-location upstream of IA
      YUA(I) = VY(IA)-LY ! y-location upstream of IA

  ! If B is upstream
      XUB(I) = VX(IB)+LX ! x-location upstream of IB
      YUB(I) = VY(IB)+LY ! y-location upstream of IB
 
  ! Finding the distance between the upstream-points and all nodes in the domain
         DO I2=1,M
          dstA(I2)=SQRT((VX(I2)-XUA(I))**2.0_SP+(VY(I2)-YUA(I))**2.0_SP)
          dstB(I2)=SQRT((VX(I2)-XUB(I))**2.0_SP+(VY(I2)-YUB(I))**2.0_SP)
         END DO

  ! Choosing the node closest to the imaginary upstream node
         minA = MINVAL(dstA)
         minB = MINVAL(dstB)

         DO I2=1,M
           if (dstA(I2).EQ.minA) then
              Anear_node(I) = I2
           end if
           if (dstB(I2).EQ.minB) then
              Bnear_node(I) = I2
           end if
         END DO

  ! Storing the distance from the upstreamnode to the nearest node
  ! ------------------------------------------------------------
  ! A-node is upstream
      XUAdist(I) = XUA(I) - VX(Anear_node(I))
      YUAdist(I) = YUA(I) - VY(Anear_node(I))

  ! B-node is upstream
      XUBdist(I) = XUB(I) - VX(Bnear_node(I))
      YUBdist(I) = YUB(I) - VY(Bnear_node(I))
  END DO 

   RETURN
   END SUBROUTINE SETUP_TVD

END MODULE MOD_TVD
