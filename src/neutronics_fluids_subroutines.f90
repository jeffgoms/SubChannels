        SUBROUTINE FISSION_MATRIX_VEC(PSI,SIGMA_F_MATRIX_KEEP,S_NEUTRONS, NX,NY,NZ,NG)
! MATRIX VECTOR FOR S_NEUTRONS=FISSION_MATRIX*PSI
        IMPLICIT NONE
        INTEGER, intent( in ) :: NX,NY,NZ,NG
        REAL, intent( in ) :: PSI(NX,NY,NZ,NG),SIGMA_F_MATRIX_KEEP(NX,NY,NZ,NG,NG)
        REAL, intent( inout ) :: S_NEUTRONS(NX,NY,NZ,NG)
! Local variables...
        INTEGER I,J,K,G,GD

        S_NEUTRONS=0.0
! MATRIX VECTOR FOR S_NEUTRONS=FISSION_MATRIX*PSI
        DO GD=1,NG ! No of energy groups
        DO G=1,NG ! No of energy groups

        DO K=2,NZ-1
        DO J=2,NY-1
        DO I=2,NX-1
! Matrix vector...
           S_NEUTRONS(I,J,K,G)=S_NEUTRONS(I,J,K,G)  + SIGMA_F_MATRIX_KEEP(I,J,K,G,GD)*PSI(I,J,K,GD) 
        END DO
        END DO
        END DO

        END DO ! DO G=1,NG
        END DO ! DO GD=1,NG
        RETURN
        END SUBROUTINE FISSION_MATRIX_VEC




        SUBROUTINE NEUTRONICS_STEP(TF,TFOLD, TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P,VEL,VELOLD, &
                             TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS, P_BCS_ON, &
                             DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NG,NR,NPHASE, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, RADIAL_ROD_MATERIAL, &
                             TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, &
                             TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, &
                             PSI,PSIOLD,C,COLD,S_NEUTRONS,SPEED,NEUTRON_DIFF,SIGMA_MATRIX, NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER,NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G, &
                             NDELAY,BETA_EFF,BETA,LAMBDA,XI_D,XI_P,SIGMA_F_MATRIX,SIGMA_F,NU_F, &
                             NSTEP_ITS, & 
                             INTERP_TEMPERATURE, NEUTRON_DIFF_TEMPERATURE,SIGMA_A_TEMPERATURE, SIGMA_S_TEMPERATURE, SIGMA_F_TEMPERATURE, &
                             XI_P_TEMPERATURE,XI_D_TEMPERATURE,NU_F_TEMPERATURE, MATERIAL,NMATERIAL, NX_T1,NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP, &
                             VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                             VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER, NTIME)
! This sub calculates the new Fuel temp TFUEL (also TFOLD) and water temperatures TW (also TWOLD) and fluids vel, p, volfra if NPHASE>1. 
! VEL is the water velocity in the vertical direction -can be +ve or -ve.
!
! Axis-symmetric diffusion of temperature for rods...
! TF,TFOLD are the lattest and previos time level values of rod temperatures. 
! TW,TWOLD are the lattest and previos time level values of water temperatures. 
! VOLFRA are the phase volume fractions.
! TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS contains the bcs of water temp, water density, velocity, volume fraction and pressure.
! P_BCS_ON =1 switches on the pressure bc for a particular phase, P_BCS_ON =0 assume that we have a vel b.c.
! DEN_CP_FUEL_RODS =CP*ROW of the fuel rods (could be control rods and includes cladding materials). 
! DIFF_FUEL_RODS = conduction for fuel/control rods (could be control rods and includes cladding materials) - is input to this sub. 
! ROD_RADIUS_NODES are the nodes of the cells (side of the cells) for the model of the fuel pin in the radial direction.
! If one need a non-uniform mesh spacing in radial direction then put this adjustment into the diffusion coefficient 
! in the radial direction. 
!
! x-coord will be along x-direction. 
! y-coord will be along radial direction.  x-direction. 
! z-coord will be along z-direction. 
! g-coord will be along radial direction. 
! NR is the no of cells in the radial direction of the fuel rod
! NX,NY,NZ no of cells in the x,y,z directions and the no of fuel rods = NX*NY (THIS INCLUDES HALOS)
! SOLVER OPTIONS:
! The smaller TF_RELAX the better the convergence =0,1,2,3 
! TF_RELAX = 0 : tridiagonal Gauss-Seidel.  
! TF_RELAX = 1 : pt Gauss-Seidel.  
! TF_RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems. 
! TF_RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) 
!   (only use for small nos of energy groups<10 as we invert a NGxNG matrix). - great for interphase coupling. 
! TF_MAX_ITS = max no of iterations.
! TF_ERROR_TOLER = error tolerence for convergence mas=x difference between temps from 2 consecutive iterations. 
! Similarly for other solver options. 
         IMPLICIT NONE
         REAL, PARAMETER :: HEAT_PER_FISSION = 3.2E-11 ! Joules per fission event.
         INTEGER, intent( in ) :: NX,NY,NZ,NG,NR,NDELAY,NPHASE, TF_MAX_ITS, TW_MAX_ITS, NSTEP_ITS
         REAL, intent( inout ) :: PSI(NX,NY,NZ,NG), C(NX,NY,NZ,NDELAY)
         REAL, intent( in ) :: PSIOLD(NX,NY,NZ,NG), COLD(NX,NY,NZ,NDELAY), S_NEUTRONS(NX,NY,NZ,NG)
         REAL, intent( in ) :: DX(3),ROD_RADIUS_NODES(NR+1),DT
         REAL, intent( inout ) :: SPEED(NG)
         REAL, intent( inout ) :: TF(NX,NY,NZ,NR),TFOLD(NX,NY,NZ,NR), TW(NX+1,NY+1,NZ,NPHASE),TWOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: VOLFRA(NX+1,NY+1,NZ,NPHASE),VOLFRAOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: DEN_WATER(NX+1,NY+1,NZ,NPHASE),DEN_WATER_OLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: P(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: VEL(NX+1,NY+1,NZ+1,NPHASE),VELOLD(NX+1,NY+1,NZ+1,NPHASE)
         REAL, intent( in ) :: TW_BCS(NX+1,NY+1,NZ,NPHASE), DEN_WATER_BCS(NX+1,NY+1,NZ,NPHASE) 
         REAL, intent( in ) :: VEL_BCS(NX+1,NY+1,NZ+1,NPHASE), VOLFRA_BCS(NX+1,NY+1,NZ,NPHASE), P_BCS(NX+1,NY+1,NZ,NPHASE)
         INTEGER, intent( in ) :: P_BCS_ON(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( in ) :: DEN_CP_FUEL_RODS(NX,NY,NZ,NR), DIFF_FUEL_RODS(NX,NY,NZ,NR)
         INTEGER, intent( in ) :: RADIAL_ROD_MATERIAL(NX,NY,NZ,NR)
         REAL, intent( in ) :: TF_ERROR_TOLER,TW_ERROR_TOLER
         INTEGER, intent( in ) :: TF_RELAX,TW_RELAX
         REAL, intent( inout ) :: NEUTRON_DIFF(NX,NY,NZ,NG)
         REAL, intent( inout ) :: SIGMA_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F_MATRIX(NX,NY,NZ,NG,NG), SIGMA_F(NX,NY,NZ,NG)
         REAL, intent( inout ) :: XI_P(NX,NY,NZ,NG),XI_D(NX,NY,NZ,NG),NU_F(NX,NY,NZ,NG)
         REAL, intent( in ) :: BETA_EFF,LAMBDA(NDELAY),BETA(NDELAY) 
         INTEGER, intent( in ) :: NX_T1, NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP
         INTEGER, intent( in ) :: MATERIAL(NX,NY,NZ),NMATERIAL
         REAL, intent( in ) :: INTERP_TEMPERATURE(MAX_NX_TI,NDIM_INTERP )
         REAL, intent( in ) :: NEUTRON_DIFF_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2)
        REAL, intent( in ) ::  SIGMA_A_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),SIGMA_S_TEMPERATURE(NG,NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL, intent( in ) :: SIGMA_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL, intent( in ) :: XI_P_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),XI_D_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL, intent( in ) :: NU_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL, intent( in ) :: VOLFRA_ERROR_TOLER, VEL_ERROR_TOLER, NEU_ERROR_TOLER, NEU_ERROR_TOLER_IN_G
         INTEGER, intent( in ) :: VOLFRA_MAX_ITS, VEL_MAX_ITS, NEU_MAX_ITS, NEU_MAX_ITS_IN_G
         INTEGER, intent( in ) :: VOLFRA_RELAX, VEL_RELAX, NEU_RELAX
         INTEGER, intent( in ) :: NTIME
! Local variables...
         INTEGER :: ITS,I,J,K,IR
         REAL :: POWER_PER_ROD, RNORM, R_WEIGHT
         REAL :: RFISS(NR)
         REAL, DIMENSION(:,:, :,:), allocatable :: S_FUEL_RODS

        DO ITS=1,NSTEP_ITS  !  ITERATE WITHIN A TIME STEP IF NEEDED. 

           IF((TW_MAX_ITS.GT.0).AND.(NTIME.NE.0)) THEN
! work out what the power source is
              ALLOCATE(S_FUEL_RODS(NX,NY,NZ,NR)); S_FUEL_RODS=0.0
              DO K=2,NZ-1
              DO J=2,NY-1
              DO I=2,NX-1
                 POWER_PER_ROD = HEAT_PER_FISSION * SUM( SIGMA_F(I,J,K,:)*PSI(I,J,K,:)) ! Watts/cm^3. 
           
                 RFISS(:) = REAL(MAX(2-RADIAL_ROD_MATERIAL(I,J,K,:),0))! Lets say that IMAT=1 is fissile fuel rods. 
                 RNORM=0.0
                 DO IR=2,NR-1
                    R_WEIGHT=RFISS(I)*( ROD_RADIUS_NODES(I+1)**2 - ROD_RADIUS_NODES(I)**2 )
                    RNORM = RNORM + R_WEIGHT
                    S_FUEL_RODS(I,J,K,IR) = R_WEIGHT*POWER_PER_ROD
                 END DO
      
                 S_FUEL_RODS(I,J,K,:) = S_FUEL_RODS(I,J,K,:)/ MAX(1.E-10,RNORM) 
              END DO
              END DO
              END DO
           CALL TEMP_CALC(TF,TFOLD, TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P,VEL,VELOLD, &
                          TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS, P_BCS_ON, &
                          DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NR,NPHASE, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, S_FUEL_RODS, &
                          TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, &
                          TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, &
                          VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                          VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER)
           ENDIF ! IF(TW_MAX_ITS.GT.0) THEN

           CALL INTERPOLATE_XSECTIONS_4D(TF,TW,VOLFRA,DEN_WATER,NEUTRON_DIFF,SIGMA_MATRIX,SIGMA_F_MATRIX,SIGMA_F, XI_P,XI_D,NU_F, ROD_RADIUS_NODES, &
                            NX,NY,NZ,NG,NR,NDELAY,NPHASE, &
                            INTERP_TEMPERATURE, NEUTRON_DIFF_TEMPERATURE,SIGMA_A_TEMPERATURE, SIGMA_S_TEMPERATURE, SIGMA_F_TEMPERATURE, &
                            XI_P_TEMPERATURE,XI_D_TEMPERATURE,NU_F_TEMPERATURE, MATERIAL,NMATERIAL, NX_T1,NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP, &
                            BETA_EFF, NTIME )

! Caculate scalar flux PSI:
           CALL NEUTRONICS_DIFFUSION(PSI,PSIOLD,S_NEUTRONS,SPEED,NEUTRON_DIFF,SIGMA_MATRIX,DT,DX, NX,NY,NZ,NG, &
                                     NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER,NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G, &
                                     NDELAY,BETA_EFF,BETA,LAMBDA,XI_D,XI_P,SIGMA_F_MATRIX,SIGMA_F,C) 

! Delayed neutron calculation for C:
           CALL CALC_DELAYED_PRECURSORS(C,COLD,PSI, NX,NY,NZ,NG,NDELAY, LAMBDA,BETA, NU_F,SIGMA_F, DT)

        END DO ! END OF DO ITS=1,NSTEP_ITS
        RETURN
        END SUBROUTINE NEUTRONICS_STEP




        SUBROUTINE INTERPOLATE_XSECTIONS_4D(TF,TW,VOLFRA,DEN_WATER,NEUTRON_DIFF,SIGMA_MATRIX,SIGMA_F_MATRIX,SIGMA_F, XI_P,XI_D,NU_F, ROD_RADIUS_NODES, &
                            NX,NY,NZ,NG,NR,NDELAY,NPHASE, &
                            INTERP_TEMPERATURE, NEUTRON_DIFF_TEMPERATURE,SIGMA_A_TEMPERATURE, SIGMA_S_TEMPERATURE, SIGMA_F_TEMPERATURE, &
                            XI_P_TEMPERATURE,XI_D_TEMPERATURE,NU_F_TEMPERATURE, MATERIAL,NMATERIAL, NX_T1,NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP, &
                            BETA_EFF, NTIME )
! This subroutine does 4D interpolations of the cross sections in say:
! ---mean radial temperature of in fuel/control rod
! ---water temperature
! ---water voidage
! ---central axies of rod temperature of fuel/control rod.  
! One could not use the last two, say, and set NX_T3=1,NX_T4=1.
! NX,NY,NZ,NG no of cells in x,y,z directions, Number of energy groups, 
! NR=no of radial cells in fuel rods. 
! NDELAY=no of delayed neutrons
! NX_T1, NX_T2, NX_T3,NX_T4 are the 4D dimensions used for temperature interpolation of x-sections. 
        IMPLICIT NONE
        INTEGER, intent( in ) :: NX,NY,NZ,NG,NR,NDELAY,NPHASE, NMATERIAL
        REAL, intent( in ) :: TF(NX,NY,NZ,NR),TW(NX+1,NY+1,NZ,NPHASE),VOLFRA(NX+1,NY+1,NZ,NPHASE),DEN_WATER(NX+1,NY+1,NZ,NPHASE)
        REAL, intent( inout ) :: NEUTRON_DIFF(NX,NY,NZ,NG)
        REAL, intent( inout ) :: SIGMA_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F(NX,NY,NZ,NG)
        REAL, intent( inout ) :: XI_P(NX,NY,NZ,NG),XI_D(NX,NY,NZ,NG),NU_F(NX,NY,NZ,NG)
        REAL, intent( in ) :: ROD_RADIUS_NODES(NR+1)
        INTEGER, intent( in ) :: NX_T1, NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP
        REAL, intent( in ) :: INTERP_TEMPERATURE(MAX_NX_TI,NDIM_INTERP )
        REAL, intent( in ) :: NEUTRON_DIFF_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
        REAL, intent( in ) ::  SIGMA_A_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),SIGMA_S_TEMPERATURE(NG,NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
        REAL, intent( in ) :: SIGMA_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
        REAL, intent( in ) :: XI_P_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),XI_D_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
        REAL, intent( in ) :: NU_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
        INTEGER, intent( in ) :: MATERIAL(NX,NY,NZ) 
        REAL, intent( in ) :: BETA_EFF
        INTEGER, intent( in ) :: NTIME
! Local variables...
        REAL, DIMENSION( :, :, :, : ), allocatable :: TWF_FOR_INTERPOLATION
        REAL, DIMENSION( : ), allocatable :: MID_RADIUS
        REAL :: ONEIN(NDIM_INTERP,2)
        INTEGER :: I,J,K,G,GD,IDIM, II,JJ,KK,MM, IR, IONLY_ONE(NDIM_INTERP), IMAT
        INTEGER :: IP,JP,KP,MP
        REAL :: R_WEIGHT, RNORM, RBETA_EFF                            
        INTEGER :: NX_TEMPERATURE(NDIM_INTERP), MAT(NDIM_INTERP)

        ALLOCATE(TWF_FOR_INTERPOLATION(NX,NY,NZ,NDIM_INTERP), MID_RADIUS(NR) )
        DO IR=1,NR
           MID_RADIUS(IR)=0.5*(ROD_RADIUS_NODES(IR)+ROD_RADIUS_NODES(IR+1))
        END DO

        RNORM = SUM(MID_RADIUS(:))
        DO K=2,NZ-1
        DO J=2,NY-1
        DO I=2,NX-1
           TWF_FOR_INTERPOLATION(I,J,K,1) = 0.25*( TW(I,J,K,1) + TW(I+1,J,K,1) + TW(I,J+1,K,1) + TW(I+1,J+1,K,1) ) ! water temp
           TWF_FOR_INTERPOLATION(I,J,K,2) = SUM( TF(I,J,K,:)*MID_RADIUS(:) )/RNORM ! mean rod temp
           IF(NPHASE.GT.1) THEN
              TWF_FOR_INTERPOLATION(I,J,K,3) = 1.0-VOLFRA(I,J,K,1) ! ASSUME THE 1ST PHASE IS WATER (THIS IS VOID FRACTION)
              TWF_FOR_INTERPOLATION(I,J,K,4) = TF(I,J,K,1) ! central axies temperature of fuel/control rod. 
           ELSE
              TWF_FOR_INTERPOLATION(I,J,K,3) = 0.0
              TWF_FOR_INTERPOLATION(I,J,K,4) = 0.0
           ENDIF
        END DO
        END DO
        END DO
        
        NX_TEMPERATURE(1)=NX_T1; NX_TEMPERATURE(2)=NX_T2; NX_TEMPERATURE(3)=NX_T3; NX_TEMPERATURE(4)=NX_T4
        IONLY_ONE(:)= 1 - MAX(0,2-NX_TEMPERATURE(:))

        DO K=2,NZ-1
        DO J=2,NY-1
        DO I=2,NX-1
           IMAT = MATERIAL(I,J,K) ! NMATERIAL is the no of different types of materials we have e.g. fuel, cladding
           DO IDIM=1,NDIM_INTERP
              CALL TEMMATINTERP(MAT(IDIM), ONEIN(IDIM,1), ONEIN(IDIM,2), &
                      TWF_FOR_INTERPOLATION(I,J,K,IDIM), NX_TEMPERATURE(IDIM), INTERP_TEMPERATURE(:,IDIM))
! NB ONEIN(IDIM,2) = 1.-ONEIN(IDIM,1)
           END DO
           II=MAT(1); JJ=MAT(2); KK=MAT(3); MM=MAT(4)

           NEUTRON_DIFF(I,J,K,:) =0.0
           SIGMA_MATRIX(I,J,K,:,:) = 0.0
           SIGMA_F(I,J,K,:) = 0.0
           XI_P(I,J,K,:) = 0.0
           XI_D(I,J,K,:) = 0.0
           NU_F(I,J,K,:) = 0.0

           DO IP=0,1*IONLY_ONE(1)
           DO JP=0,1*IONLY_ONE(2)
           DO KP=0,1*IONLY_ONE(3)
           DO MP=0,1*IONLY_ONE(4)

              R_WEIGHT= ONEIN(1,IP+1) * ONEIN(2,JP+1) * ONEIN(3,KP+1) * ONEIN(4,MP+1)

              NEUTRON_DIFF(I,J,K,:)= NEUTRON_DIFF(I,J,K,:)  +  R_WEIGHT * NEUTRON_DIFF_TEMPERATURE(:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP) 
              SIGMA_MATRIX(I,J,K,:,:) = SIGMA_MATRIX(I,J,K,:,:) + R_WEIGHT * SIGMA_S_TEMPERATURE(:,:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP)
              DO G=1,NG
                 SIGMA_MATRIX(I,J,K,G,G) = SIGMA_MATRIX(I,J,K,G,G) + R_WEIGHT * SIGMA_A_TEMPERATURE(G,IMAT,II+IP,JJ+JP,KK+KP,MM+MP)
              END DO
              SIGMA_F(I,J,K,:) = SIGMA_F(I,J,K,:) + R_WEIGHT * SIGMA_F_TEMPERATURE(:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP) 

! Prompt fission spectrum...
              XI_P(I,J,K,:)= XI_P(I,J,K,:) + R_WEIGHT * XI_P_TEMPERATURE(:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP) 
! Delayed fission spectrum...
              XI_D(I,J,K,:)= XI_D(I,J,K,:) + R_WEIGHT * XI_D_TEMPERATURE(:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP) 
! No of neutrons produced per fission event for each group...
              NU_F(I,J,K,:)= NU_F(I,J,K,:) +  R_WEIGHT * NU_F_TEMPERATURE(:,IMAT,II+IP,JJ+JP,KK+KP,MM+MP) 

           END DO
           END DO
           END DO
           END DO

        END DO
        END DO
        END DO ! DO K=2,NZ-1

        RBETA_EFF=0.0
        IF(NTIME.GE.1) RBETA_EFF=BETA_EFF
        DO GD=1,NG
        DO G=1,NG
           SIGMA_F_MATRIX(:,:,:,G,GD) =  (1.0-RBETA_EFF)*XI_P(:,:,:,G)*NU_F(:,:,:,GD)*SIGMA_F(:,:,:,GD) 
        END DO
        END DO

        RETURN
        END SUBROUTINE INTERPOLATE_XSECTIONS_4D




    SUBROUTINE TEMMATINTERP(MAT, ONEIN, TWOIN, &
         TEMP, NDSIDI1, TEMMAT)
      ! Interpolate temps to find MAT and 
      ! interpolation parameters ONEIN,TWOIN...
      !use FLDebug
      IMPLICIT NONE
      INTEGER, intent(inout) :: MAT
      REAL, intent(inout) :: ONEIN, TWOIN
      REAL, intent(in) :: TEMP
      INTEGER, intent(in) :: NDSIDI1
      REAL, DIMENSION(NDSIDI1), intent(in) :: TEMMAT
      ! Local
      INTEGER :: JMAT

      IF( NDSIDI1.EQ.1) THEN ! only one dimension in this space...

         MAT=1
         ONEIN = 1.
         TWOIN = 0.
      
      ELSE

         IF( TEMP >= TEMMAT(NDSIDI1) ) THEN
            MAT = NDSIDI1 - 1
            ONEIN = 0.
            TWOIN = 1.

         ELSE IF( TEMP <= TEMMAT(1) ) THEN 
            MAT = 1
            ONEIN = 1.
            TWOIN = 0.

         ELSE

            DO JMAT = NDSIDI1 - 1, 1, -1

               IF( (TEMP >= TEMMAT(JMAT)) .AND. &
                   (TEMP <= TEMMAT(JMAT + 1)) ) MAT = JMAT

            END DO

            TWOIN =  (TEMP - TEMMAT(MAT)) / (TEMMAT(MAT + 1) - TEMMAT(MAT))
            ONEIN = 1. - TWOIN

         ENDIF

      ENDIF ! IF( NDSIDI1.EQ.1) THEN 

      RETURN
    END SUBROUTINE TEMMATINTERP




        SUBROUTINE CALC_DELAYED_PRECURSORS(C,COLD,PSI, NX,NY,NZ,NG,NDELAY, LAMBDA,BETA, NU_F,SIGMA_F, DT)
! Calculate the delayed neutrons C. 
! C on the same grid as PSI - the flux
        IMPLICIT NONE
        INTEGER, intent( in ) :: NX,NY,NZ,NG,NDELAY
        REAL, intent( inout ) :: C(NX,NY,NZ,NDELAY)
        REAL, intent( in ) :: COLD(NX,NY,NZ,NDELAY)
        REAL, intent( in ) :: PSI(NX,NY,NZ,NG), SIGMA_F(NX,NY,NZ,NG)
        REAL, intent( in ) :: LAMBDA(NDELAY),BETA(NDELAY), NU_F(NX,NY,NZ,NG), DT
! Local variables...
        REAL :: DIA(NDELAY)
        INTEGER :: D,G
        
        DIA(:)=1./DT  + LAMBDA(:)

        C = COLD
        DO D=1,NDELAY
           DO G=1,NG
              C(:,:,:,D) = C(:,:,:,D) + BETA(D)*NU_F(:,:,:,G)*SIGMA_F(:,:,:,G)*PSI(:,:,:,G)/DIA(D) 
           END DO
        END DO
        RETURN
        END SUBROUTINE CALC_DELAYED_PRECURSORS




        SUBROUTINE TEMP_CALC(TF,TFOLD, TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P,VEL,VELOLD, &
                             TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS, P_BCS_ON, &
                             DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NR,NPHASE, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, S_FUEL_RODS, &
                             TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, &
                             TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, &
                             VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                             VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER)
! This sub calculates the new Fuel temp TFUEL (also TFOLD) and water temperatures TW (also TWOLD). 
! VELOCITY_UP is the water velocity in the vertical direction -can be +ve or -ve.
! Axis-symmetric diffusion of temperature for rods...
! TF,TFOLD are the lattest and previos time level values of rod temperatures. 
! TW,TWOLD are the lattest and previos time level values of water temperatures. 
! DEN_CP_FUEL_RODS =CP*ROW of the fuel rods (could be control rods and includes cladding materials) - is input to this sub. 
! DIFF_FUEL_RODS = conduction for fuel/control rods (could be control rods and includes cladding materials) - is input to this sub. 
! DEN_CP_WATER= density * CP for water and for each phase of it (water, vapour, droplets) - is input for NPHASE=1 else caculated here. 
! ROD_RADIUS_NODES are the nodes of the cells (side of the cells) for the model of the fuel pin in the radial direction.
! If one need a non-uniform mesh spacing in radial direction then put this adjustment into the diffusion coefficient 
! in the radial direction. 
! TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS contains the bcs of water temp, water density, velocity, volume fraction and pressure.
! P_BCS_ON =1 switches on the pressure bc for a particular phase, P_BCS_ON =0 assume that we have a vel b.c.
!
! x-coord will be along x-direction. 
! y-coord will be along radial direction.   
! z-coord will be along z-direction. 
! g-coord will be along radial direction. 
! NR is the no of cells in the radial direction of the fuel rod
! NX,NY,NZ no of cells in the x,y,z directions and the no of fuel rods = NX*NY (THIS INCLUDES HALOS)
! SOLVER OPTIONS:
! The smaller TF_RELAX the better the convergence =0,1,2,3 
! TF_RELAX = 0 : tridiagonal Gauss-Seidel.  
! TF_RELAX = 1 : pt Gauss-Seidel.  
! TF_RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems. 
! TF_RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) 
!   (only use for small nos of energy groups<10 as we invert a NGxNG matrix). - great for interphase coupling. 
! TF_MAX_ITS = max no of iterations.
! TF_ERROR_TOLER = error tolerence for convergence mas=x difference between temps from 2 consecutive iterations. 
! Similarly for other solver options.
         IMPLICIT NONE
         INTEGER, PARAMETER :: NITS_NONLIN_TW = 0 
! NITS_NONLIN_TW controls no of non-linear iterations of flux limiter. If NITS_NONLIN_TW = 0 have no flux limiting. 
         INTEGER, intent( in ) :: NX,NY,NZ,NR,NPHASE, TF_MAX_ITS, TW_MAX_ITS
         REAL, intent( in ) :: DX(3),ROD_RADIUS_NODES(NR+1),DT
         REAL, intent( inout ) :: TF(NX,NY,NZ,NR),TFOLD(NX,NY,NZ,NR), TW(NX+1,NY+1,NZ,NPHASE),TWOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: VOLFRA(NX+1,NY+1,NZ,NPHASE),VOLFRAOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: DEN_WATER(NX+1,NY+1,NZ,NPHASE),DEN_WATER_OLD(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
         REAL, intent( inout ) :: VEL(NX+1,NY+1,NZ+1,NPHASE),VELOLD(NX+1,NY+1,NZ+1,NPHASE)
         REAL, intent( in ) :: TW_BCS(NX+1,NY+1,NZ,NPHASE), DEN_WATER_BCS(NX+1,NY+1,NZ,NPHASE) 
         REAL, intent( in ) :: VEL_BCS(NX+1,NY+1,NZ+1,NPHASE), VOLFRA_BCS(NX+1,NY+1,NZ,NPHASE), P_BCS(NX+1,NY+1,NZ,NPHASE)
         INTEGER, intent( in ) :: P_BCS_ON(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( in ) :: DEN_CP_FUEL_RODS(NX,NY,NZ,NR), DIFF_FUEL_RODS(NX,NY,NZ,NR), S_FUEL_RODS(NX,NY,NZ,NR)
         REAL, intent( in ) :: TF_ERROR_TOLER,TW_ERROR_TOLER
         INTEGER, intent( in ) :: TF_RELAX,TW_RELAX
         REAL, intent( in ) :: VOLFRA_ERROR_TOLER, VEL_ERROR_TOLER
         INTEGER, intent( in ) :: VOLFRA_MAX_ITS, VEL_MAX_ITS
         INTEGER, intent( in ) :: VOLFRA_RELAX, VEL_RELAX
! Local variables...
         REAL :: DX_TEMP(3), INCOME
         INTEGER :: NX_TEMP,NY_TEMP,NZ_TEMP,NR_TEMP
         REAL, DIMENSION( :,:, :,: ), allocatable :: TW2, CP_WATER, VOL_DEN_CP_WATER
         REAL, DIMENSION( :,:, :,: ), allocatable :: SIGMA_W, SIGMA_ROD
         REAL, DIMENSION( :,:, :,: ), allocatable :: SW
         REAL, DIMENSION( :,:, :,:, : ), allocatable :: HW
         REAL, DIMENSION(:,:,:,:,:,:), allocatable :: KAPPAW
         REAL, DIMENSION(:,:, :,:, :), allocatable :: CONSERV_VERT_ADV
         INTEGER :: I,J,K,G, IUP, IPHASE, ITS_NONLIN_TW

! Calculate CP_WATER for the different water phases...
! *************************************************
        print *, 'Calculate CP'
        ALLOCATE(CP_WATER(NX+1,NY+1,NZ,NPHASE), VOL_DEN_CP_WATER(NX+1,NY+1,NZ,NPHASE))
        CALL CALCULATE_CP_WATER(DEN_WATER, CP_WATER, TW, P, NX,NY,NZ,NPHASE)
        VOL_DEN_CP_WATER=VOLFRA*DEN_WATER*CP_WATER

! Calculate SIGMA_W & SIGMA_ROD for both water and fuel rods...
! *************************************************
        ALLOCATE(SIGMA_W(NX+1,NY+1,NZ,NPHASE), SIGMA_ROD(NX+1,NY+1,NZ,NPHASE)) 

        SIGMA_W = 0.0; SIGMA_ROD = 0.0

        print *, 'calling CALC_SIGMA_ROD_WATER'  
        CALL CALC_SIGMA_ROD_WATER( SIGMA_W, SIGMA_ROD, TW, TF, VEL, DX,ROD_RADIUS_NODES, NX,NY,NZ,NPHASE,NR )
        print *, 'return from CALC_SIGMA_ROD_WATER'

! **********************************************************************************
! Apply bcs to TW and other fields at top and bottom ******************************* 
        print *, 'apply bcs'
        K=1
        TW(:,:,K,:)=TW_BCS(:,:,K,:) ! This is only for advection and active only if vel comes into domain from below. 
        DEN_WATER(:,:,K,:) = DEN_WATER_BCS(:,:,K,:)
        VOLFRA(:,:,K,:) = VOLFRA_BCS(:,:,K,:)
        P(:,:,K) = P_BCS(:,:,K,1)
        VEL(:,:,K,:) = VEL_BCS(:,:,K,:) ! need this for momentum comming into domain even if pressure bc
        K=NZ
        TW(:,:,K,:)=TW_BCS(:,:,K,:) ! This is only for advection and active only if vel comes into domain from top.
        DEN_WATER(:,:,K,:) = DEN_WATER_BCS(:,:,K,:)
        VOLFRA(:,:,K,:) = VOLFRA_BCS(:,:,K,:)
        P(:,:,K) = P_BCS(:,:,K,1)
        VEL(:,:,K+1,:) = VEL_BCS(:,:,K+1,:) ! need this for momentum comming into domain even if pressure bc

! Axis-symmetric diffusion calculation for rods...
! *************************************************
        CALL AXIS_SYM_DIFF_FUEL_RODS(TF,TFOLD, TW,TWOLD, SIGMA_ROD, DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NR,NPHASE, &
                                     S_FUEL_RODS, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER)


! Calculate the water temperature...
! *************************************************
        IF(NPHASE.GT.1) THEN ! multi-phase fluids calculation FOR VELOCITY, PRESSURE, VOLUME FRACTION, DEN_WATER...
           CALL FLUIDS_CALC(TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P, VEL,VELOLD, &
                            VEL_BCS, P_BCS_ON, &
                            ROD_RADIUS_NODES,DX,DT, NX,NY,NZ,NR,NPHASE, &
                            VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                            VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER)
        ENDIF

! calculate HW,SW,KAPPAW for the water subchannels rods...
        ALLOCATE(SW(NX+1,NY+1,NZ,NPHASE))
        ALLOCATE(HW(NX+1,NY+1,NZ,NPHASE,NPHASE))
        DO IPHASE=1,NPHASE ! No of fluid phases
           SW(:,:,:,IPHASE)=VOL_DEN_CP_WATER(:,:,:,IPHASE)*TWOLD(:,:,:,IPHASE)/DT  + SIGMA_W(:,:,:,IPHASE)*TF(:,:,:,NR-1) 
           HW(:,:,:,IPHASE,IPHASE) = VOL_DEN_CP_WATER(:,:,:,IPHASE)/DT  + SIGMA_W(:,:,:,IPHASE)
        END DO
        print *, "max", maxval(HW)
        print *, "HW", HW
        print *, "sigma_W", SIGMA_W

! Add in the vertical velocity (In non-conservative form) as a diffusion term...
        ALLOCATE(KAPPAW(2,3,NX+1,NY+1,NZ,NPHASE)); KAPPAW=0.0 ! set diffusion to zero until we have a model for this.

        ALLOCATE(CONSERV_VERT_ADV(2,NX+1,NY+1,NZ,NPHASE),TW2(NX+1,NY+1,NZ,NPHASE))
        DO ITS_NONLIN_TW = 1, MAX(1,NITS_NONLIN_TW)
           CONSERV_VERT_ADV=0.0
           DO K=2,NZ-1
           DO IUP=1,2 
              CONSERV_VERT_ADV(IUP,:,:,K,:) = VOL_DEN_CP_WATER(:,:,K,:) * VEL(:,:,K-1+IUP,:) /DX(3) 
           END DO ! DO IUP=1,2
           END DO
! apply limiting to TW
           IF(NITS_NONLIN_TW.NE.0) CALL AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV, VEL, TW, NX+1,NY+1,NZ,NPHASE)
        
! CALCULATE THE TEMP OF THE WATER
        print *, "TW_RELAX", TW_RELAX
        print *, 'Calculate water temperature - diffusion_adv_cal'
           CALL DIFFUSION_ADV_CAL(TW,TW2,KAPPAW,CONSERV_VERT_ADV,CONSERV_VERT_ADV,HW,SW, &
                         NX+1,NY+1,NZ,NPHASE, .FALSE., TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, 1, 0.0)
        print *, "End of TEMP_CALC"
        END DO ! DO ITS_NONLIN_TW = 1, MAX(1,NITS_NONLIN_TW)
        RETURN
        END SUBROUTINE TEMP_CALC




        SUBROUTINE AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV, VEL, T, NX,NY,NZ,NPHASE)
! Amend CONSERV_VERT_ADV to include flux limiting of T. Can repeadely apply this for density and T say. 
        IMPLICIT NONE
        INTEGER, intent( in ) :: NX,NY,NZ,NPHASE
        REAL, DIMENSION( 2,NX,NY,NZ,NPHASE ), intent( inout ) :: CONSERV_VERT_ADV
        REAL, DIMENSION( NX,NY,NZ+1,NPHASE ), intent( in ) :: VEL
        REAL, DIMENSION( NX,NY,NZ,NPHASE ), intent( in ) :: T
! local variables...
        REAL, DIMENSION( :,:, :,: ), allocatable :: T_RATIO
        INTEGER :: K
! 
        ALLOCATE(T_RATIO(NX,NY,NZ+1,NPHASE))
! Flux limiting...
        CALL ONE_D_LIMIT(T_RATIO, VEL, T, NX,NY,NZ,NPHASE) ! Calculate T_RATIO used for limiting.
! No flux limiting T_RATIO=1.0

        DO K=2,NZ-1 
           CONSERV_VERT_ADV(1,:,:,K,:) = CONSERV_VERT_ADV(1,:,:,K,:)*T_RATIO(:,:,K,:) 
           CONSERV_VERT_ADV(2,:,:,K,:) = CONSERV_VERT_ADV(2,:,:,K,:)*T_RATIO(:,:,K+1,:) 
        END DO

        RETURN
        END SUBROUTINE AMEND_FOR_LIMITING_CONSERV_VERT_ADV
        



        SUBROUTINE ONE_D_LIMIT(T_RATIO, VEL, T, NX,NY,NZ,NG) 
    ! This sub calculates the limited face values T_LIMIT of variable T. 
    ! VEL is an input and contains the direction of velocity (+ve or -ve).  
    ! Also returned is the ratio T_RATIO which is: the limited value/ the upwind value (=1 for upwinding).  
    ! Only return the ratio for the time being. 
    !--------------------------------------------------- 
    implicit none
    INTEGER, intent( in ) :: NX,NY,NZ,NG
    REAL, DIMENSION( NX,NY,NZ+1,NG ), intent( inout ) :: T_RATIO
    REAL, DIMENSION( NX,NY,NZ+1,NG ), intent( in ) :: VEL
    REAL, DIMENSION( NX,NY,NZ,NG ), intent( in ) :: T
    ! Local variables
    REAL, PARAMETER :: TOLER=1.0E-10, XI_LIMIT = 2.0 ! defines TVD curve on the NVD diagram. 
    REAL :: DENOIN,CTILIN,INCOME, DENOOU,CTILOU,FTILIN,FTILOU, T_LIMIT
    INTEGER :: I,J,K,G
! cell numbering for +ve vel:
! k-1 ! k !  k+1
!   U ! C !f D
! \hat I=(I -I_U)/(I_D - I_U)  

       ! Calculate near the boundaries
!      DO G=1,NG
!      DO J=2,NY-1
!      DO I=2,NX-1
!         K=1
!         INCOME = 0.5*(SIGN(1.0, VEL(I,J,K+1,G) ) +1.0)
!         T_LIMIT(I,J,K+1,G) =         INCOME * T(I,J,K,G) + (1.0 - INCOME) * T(I,J,K+1,G) 
!         K=NZ-1
!         INCOME = 0.5*(SIGN(1.0, VEL(I,J,K+1,G) ) +1.0)
!         T_LIMIT(I,J,K+1,G) =         INCOME * T(I,J,K,G) + (1.0 - INCOME) * T(I,J,K+1,G) 
!      END DO
!      END DO
!      END DO ! DO G=1,NG

      
! Calculate away from the boundaries (top and bottom) 
      T_RATIO = 1.0 ! This assumes an upwind approximation. 
      DO G=1,NG
      DO K=2,NZ-2
      DO J=2,NY-1
      DO I=2,NX-1

! advection upwards
          DENOIN = SIGN( MAX(ABS( T(I,J,K+1,G) - T(I,J,K-1,G) ), TOLER), T(I,J,K+1,G) - T(I,J,K-1,G) )
          CTILIN = ( T(I,J,K,G) - T(I,J,K-1,G) ) / DENOIN
          INCOME = 0.5*(SIGN(1.0, VEL(I,J,K,G) ) +1.0)

! Calculate normalisation parameters for out going velocities
! advection downwards
! cell numbering for -ve vel:
! k-1 ! k !  k+1 ! k+2
!     ! D !f C   ! U
! \hat I=(I -I_U)/(I_D - I_U)  
          DENOOU = SIGN( MAX(ABS( T(I,J,K,G) - T(I,J,K+2,G) ), TOLER), T(I,J,K,G) - T(I,J,K+2,G) )
          CTILOU = ( T(I,J,K+1,G) - T(I,J,K+2,G) ) / DENOOU

          FTILIN = ( 0.5*(T(I,J,K,G)+T(I,J,K+1,G)) - T(I,J,K-1,G) ) / DENOIN ! the mean with 0.5 coeff is the high order value of T. 
          FTILOU = ( 0.5*(T(I,J,K,G)+T(I,J,K+1,G)) - T(I,J,K+2,G) ) / DENOOU

! Velocity is going out of element
          T_LIMIT =         INCOME   * ( T(I,J,K-1,G) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * DENOIN ) &
                  + ( 1.0 - INCOME ) * ( T(I,J,K+2,G) + MAX(  MIN(FTILOU, XI_LIMIT*CTILOU, 1.0), CTILOU) * DENOOU )
          T_RATIO(I,J,K+1,G) =  T_LIMIT / SIGN( MAX(ABS( INCOME *T(I,J,K,G)+(1.0-INCOME)*T(I,J,K+1,G) ), TOLER), INCOME *T(I,J,K,G)+(1.0-INCOME)*T(I,J,K+1,G) )

      END DO
      END DO
      END DO
      END DO ! DO G=1,NG

    RETURN

  END SUBROUTINE ONE_D_LIMIT




        SUBROUTINE CALCULATE_CP_WATER(DEN_WATER, CP_WATER, TW, P, NX,NY,NZ,NPHASE)
! Calculate the heat capacity of water and steam...
        INTEGER, intent( in ) :: NX,NY,NZ,NPHASE
        REAL, intent( in ) :: TW(NX+1,NY+1,NZ,NPHASE),DEN_WATER(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
        REAL, intent( inout ) :: CP_WATER(NX+1,NY+1,NZ,NPHASE)
        
        IF(NPHASE==1) THEN ! Just water
            CP_WATER=4.2 ! Units in cm and g.
        ELSE
        ENDIF

        RETURN
        END SUBROUTINE CALCULATE_CP_WATER




        SUBROUTINE CALC_SIGMA_ROD_WATER( SIGMA_W, SIGMA_ROD, TW, TF, VEL, DX,ROD_RADIUS_NODES, NX,NY,NZ,NPHASE,NR )
        IMPLICIT NONE
        real, parameter :: PI=4.0*atan(1.0)
        INTEGER, intent( in ) :: NX,NY,NZ,NR,NPHASE
        REAL, intent( in ) :: TW(NX+1,NY+1,NZ,NPHASE), TF(NX,NY,NZ,NR), VEL(NX+1,NY+1,NZ+1,NPHASE), DX(3), ROD_RADIUS_NODES(NR+1)
        REAL, intent( inout ) :: SIGMA_W(NX+1,NY+1,NZ,NPHASE), SIGMA_ROD(NX,NY,NZ,NPHASE)
! Local variables...
        REAL, DIMENSION( : , : , :, : ), allocatable :: VEL_ROD_MEAN
        REAL :: RADIUS,RE,NU_D,H_WR
        INTEGER :: IPHASE, II, JJ, iError, i, j, k
        REAL :: KINEMATIC_VISC_WATER, k_w, Pr
        
        allocate(VEL_ROD_MEAN(NX+1,NY+1,NZ,NPHASE)); VEL_ROD_MEAN = 0.0

        RADIUS=ROD_RADIUS_NODES(NR) ! Radius of control rod.

        DO K=2,NZ-1 ! Calc mean rod velocity...
           VEL_ROD_MEAN(:,:,K,:) = 0.5*(VEL(:,:,K,:)+VEL(:,:,K+1,:))
        END DO
        IF(NPHASE==1) THEN ! Just water - Units in cm and g
           Pr = 7.0 ! Pradle no of water
           k_w =  0.6 ! conductivity of water (W/(mK)) thus 0.6 W/(mK) becomes 0.6*0.01 W/(cm K)
!           KINEMATIC_VISC_WATER = 1.E-6 ! Water kinematic viscocity m^2/s
           KINEMATIC_VISC_WATER = 1.E-2 ! Water kinematic viscocity cm^2/s
           SIGMA_W = 0.0
           SIGMA_ROD = 0.0
           RADIUS = ROD_RADIUS_NODES(NR) 
           DO IPHASE=1,NPHASE
           DO K=2,NZ-1
           DO J=2,NY-1
           DO I=2,NX-1
              Re= VEL_ROD_MEAN(I,J,K,IPHASE)*0.01 * RADIUS*0.01/  KINEMATIC_VISC_WATER
              NU_D = 0.332* Re**0.5 *Pr**0.3333333 ! Nusselt no. 
              H_WR= ((k_w*0.01)/(RADIUS*0.01))*NU_D ! heat transfer coefficient. 
              SIGMA_W(I,J,K,IPHASE) = 2.0*PI*RADIUS*H_WR/(DX(1)*DX(2)) ! Volumetric heat tranfer coeff - adjusted for discretization of water temp. 
!              SIGMA_W(I,J,K,IPHASE) = 100000.0!*2.0*PI*RADIUS*H_WR/(DX(1)*DX(2))
              !print *, "sigma_w", 2.0*PI*RADIUS*H_WR/(DX(1)*DX(2)) 
! Share between 4 control rods the volumetric heat tranfer coeff - adjusted for discretization of rod temp. 
              DO II=0,1
              DO JJ=0,1
                 SIGMA_ROD(I+II,J+JJ,K,IPHASE) = SIGMA_ROD(I+II,J+JJ,K,IPHASE) + 0.25 * RADIUS * H_WR ! Just 1/4 of the fuel rod
              END DO
              END DO
           END DO
           END DO
           END DO
           END DO ! DO IPHASE=1,NPHASE
        ELSE ! Multi-phase system
        ENDIF

        deallocate(VEL_ROD_MEAN)
        
        RETURN
        END SUBROUTINE CALC_SIGMA_ROD_WATER




        SUBROUTINE FLUIDS_CALC(TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P, VEL,VELOLD, &
                             VEL_BCS, P_BCS_ON, &
                             ROD_RADIUS_NODES,DX,DT, NX,NY,NZ,NR,NPHASE, &
                             VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                             VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER)
! This sub calculates the new the fluids velocity VEL, pressure P and volume fraction VOLFRA. 
! VELOCITY_UP is the water velocity in the vertical direction -can be +ve or -ve.
! TW,TWOLD are the lattest and previos time level values of water temperatures. 
! ROD_RADIUS_NODES are the nodes of the cells (side of the cells) for the model of the fuel pin in the radial direction.
! If one need a non-uniform mesh spacing in radial direction then put this adjustment into the diffusion coefficient 
! in the radial direction. 
! VEL_BCS contains the bcs of water velocity (the other b.cs are impliec by specifying just before this subroutine).
! P_BCS_ON =1 switches on the pressure bc for a particular phase, P_BCS_ON =0 assume that we have a vel b.c.
!
! x-coord will be along x-direction. 
! y-coord will be along radial direction.  x-direction. 
! z-coord will be along z-direction. 
! g-coord will be along radial direction. 
! NR is the no of cells in the radial direction of the fuel rod
! NX,NY,NZ no of cells in the x,y,z directions and the no of fuel rods = NX*NY (THIS INCLUDES HALOS)
! SOLVER OPTIONS:
! The smaller VEL_RELAX the better the convergence =0,1,2,3 
! VEL_RELAX = 0 : tridiagonal Gauss-Seidel.  
! VEL_RELAX = 1 : pt Gauss-Seidel.  
! VEL_RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems. 
! VEL_RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) 
!   (only use for small nos of energy groups<10 as we invert a NGxNG matrix). - great for interphase coupling. 
! VEL_MAX_ITS = max no of iterations.
! VEL_ERROR_TOLER = error tolerence for convergence mas=x difference between temps from 2 consecutive iterations. 
! Similarly for other solver options VOLFRA_RELAX. 
         IMPLICIT NONE
         REAL, PARAMETER :: GRAVTY = 9800.0 ! In units of cm,g 
         INTEGER, PARAMETER :: NITS_NONLIN_VOLFRA = 0 ! =0 no limiting for volume fraction and otherwise = no of non-linear limiting iterations
         INTEGER, PARAMETER :: NITS_NONLIN_VEL = 0 ! =0 no limiting for velocity and otherwise = no of non-linear limiting iterations
         REAL, PARAMETER :: INFINY = 1.E+14 ! Used in the application of bcs using big sping. 
         INTEGER,  PARAMETER :: NCOLOR = 3 ! No of colours to colout the tridonal pressure matrix.
         INTEGER, intent( in ) :: NX,NY,NZ,NR,NPHASE, VOLFRA_MAX_ITS, VEL_MAX_ITS
         REAL, intent( in ) :: DX(3),ROD_RADIUS_NODES(NR+1),DT ! Rod radius needed for drag correlation. 
         REAL, intent( in ) :: TW(NX+1,NY+1,NZ,NPHASE),TWOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( in ) :: VOLFRA_ERROR_TOLER, VEL_ERROR_TOLER
         INTEGER, intent( in ) :: VOLFRA_RELAX, VEL_RELAX
         REAL, intent( inout ) :: VOLFRA(NX+1,NY+1,NZ,NPHASE),VOLFRAOLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: DEN_WATER(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
         REAL, intent( in ) :: DEN_WATER_OLD(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: VEL(NX+1,NY+1,NZ+1,NPHASE),VELOLD(NX+1,NY+1,NZ+1,NPHASE)
         REAL, intent( in ) :: VEL_BCS(NX+1,NY+1,NZ+1,NPHASE)
         INTEGER, intent( in ) :: P_BCS_ON(NX+1,NY+1,NZ,NPHASE)
! Local variables...
         REAL :: DX_TEMP(3)
         INTEGER :: NX_TEMP,NY_TEMP,NZ_TEMP,NR_TEMP
         REAL, DIMENSION( : , : , :, : ), allocatable :: VOLFRA2, VEL2, CT_VEL, VEC_COL_LONG
         REAL, DIMENSION( : , : , :, : ), allocatable :: VEC_COL
         REAL, DIMENSION( : , : , : ), allocatable :: DP, PCOLOR
         REAL, DIMENSION(:,:, :,:, :,:), allocatable :: KAPPA_VOLFRA, KAPPA_VEL
         REAL, DIMENSION(:,:, :,:, :), allocatable :: CONSERV_VERT_ADV_VOLFRA, CONSERV_VERT_ADV_VOLFRA_SIGN
         REAL, DIMENSION(:,:, :,:, :), allocatable :: CONSERV_VERT_ADV_VEL, H_VOLFRA, H_VEL
         REAL, DIMENSION(:,:, :,:), allocatable :: D_DEN_D_P, DEN_WATER_NOD, VOLFRA_NOD, S_VOLFRA, S_VEL, ABC
         REAL :: MAT(NPHASE,NPHASE)
         INTEGER :: I,J,K,G,IPHASE,JPHASE,IUP,ICOLOR,JCOLOR, KK
         INTEGER :: ITS_NONLIN_VOLFRA, ITS_NONLIN_VEL


        ALLOCATE(D_DEN_D_P(NX+1,NY+1,NZ,NPHASE)) 
        CALL DEN_FUNCTION_WATER(DEN_WATER,D_DEN_D_P, TW,P, NX,NY,NZ,NPHASE) ! Use steam tables from FETCH
! Calculate density between the CVs: 
        ALLOCATE(DEN_WATER_NOD(NX+1,NY+1,NZ+1,NPHASE), VOLFRA_NOD(NX+1,NY+1,NZ+1,NPHASE)) 
        DEN_WATER_NOD(:,:,2:NZ,:) = 0.5* (DEN_WATER(:,:,1:NZ-1,:) + DEN_WATER(:,:,2:NZ,:))
        VOLFRA_NOD(:,:,2:NZ,:) = 0.5* (VOLFRA(:,:,1:NZ-1,:) + VOLFRA(:,:,2:NZ,:))

        ALLOCATE(VOLFRA2(NX+1,NY+1,NZ,NPHASE)) 
        ALLOCATE(VEL2(NX+1,NY+1,NZ+1,NPHASE)) 


! *******************************
! Solve for volume fraction...
        ALLOCATE(S_VOLFRA(NX+1,NY+1,NZ,NPHASE), H_VOLFRA(NX+1,NY+1,NZ,NPHASE,NPHASE))
        DO K=1,NZ-1
           S_VOLFRA(:,:,K,:)=DEN_WATER_OLD(:,:,K,:)*VOLFRAOLD(:,:,K,:)/DT 
        END DO
        CALL CALC_VOLFRA_MASS_EXCHANGE(H_VOLFRA, DEN_WATER, VOLFRA, VEL, P, NX,NY,NZ,NPHASE) ! Mass exchange term between the phases.
        DO IPHASE=1,NPHASE
           H_VOLFRA(:,:,:,IPHASE,IPHASE)=  H_VOLFRA(:,:,:,IPHASE,IPHASE) + DEN_WATER(:,:,:,IPHASE)/DT
        END DO
        ALLOCATE(KAPPA_VOLFRA(2,3,NX+1,NY+1,NZ,NPHASE), CONSERV_VERT_ADV_VOLFRA(2,NX+1,NY+1,NZ,NPHASE) )
        KAPPA_VOLFRA=0.0 ! Zero momentum diffusion to start with - viscocity
! the advection velocity matrix term..
        DO ITS_NONLIN_VOLFRA=1,MAX(1,NITS_NONLIN_VOLFRA) ! Non-linear iteration for limiting
           DO K=2,NZ-1
           DO IUP=1,2
              CONSERV_VERT_ADV_VOLFRA(IUP,:,:,K,:) = DEN_WATER(:,:,K-2+IUP,:)*MAX(VEL(:,:,K-1+IUP,:),0.0)/DX(3)  &
                                                   + DEN_WATER(:,:,K-1+IUP,:)*MIN(VEL(:,:,K-1+IUP,:),0.0)/DX(3) 
           END DO
           END DO
! apply limiting to convolution DEN_WATER*VOLFRA
           IF(NITS_NONLIN_VOLFRA.NE.0) CALL AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV_VOLFRA, VEL, DEN_WATER*VOLFRA, NX+1,NY+1,NZ,NPHASE)
! solve
           CALL DIFFUSION_ADV_CAL(VOLFRA,VOLFRA2,KAPPA_VOLFRA,CONSERV_VERT_ADV_VOLFRA,CONSERV_VERT_ADV_VOLFRA,H_VOLFRA,S_VOLFRA, &
                        NX+1,NY+1,NZ,NPHASE, .FALSE., VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, 1, 0.0)
        END DO ! DO ITS_NONLIN_VOLFRA=1,MAX(1,NITS_NONLIN_VOLFRA)


! *******************************
! Calculate velcity...
        ALLOCATE( S_VEL(NX+1,NY+1,NZ+1,NPHASE), H_VEL(NX+1,NY+1,NZ+1,NPHASE,NPHASE) )
        DO IPHASE=1,NPHASE
        DO K=2,NZ
           S_VEL(:,:,K,IPHASE)=DEN_WATER_NOD(:,:,K,IPHASE)*VELOLD(:,:,K,IPHASE)/DT  - (P(:,:,K)-P(:,:,K-1))/DX(3) & ! need pressure in
                         -DEN_WATER_NOD(:,:,K,IPHASE)*GRAVTY
        END DO
        END DO
        CALL CALC_VEL_ABSORPTION(H_VEL, DEN_WATER_NOD, VOLFRA_NOD, VEL, P, NX,NY,NZ,NPHASE) ! Calculate friction with pipe and between phases.
        DO IPHASE=1,NPHASE
           H_VEL(:,:,:,IPHASE,IPHASE)=  H_VEL(:,:,:,IPHASE,IPHASE) + DEN_WATER_NOD(:,:,:,IPHASE)/DT
        END DO
        ALLOCATE(KAPPA_VEL(2,3,NX+1,NY+1,NZ+1,NPHASE), CONSERV_VERT_ADV_VEL(2,NX+1,NY+1,NZ+1,NPHASE) )
        KAPPA_VEL=0.0 ! Zero momentum diffusion to start with - viscocity
! Apply Dirichlet bcs to velocity but only if not a pressure condition P_BCS_ON=0 using the big spring method. 
! Apply bcs to TW and other fields at top and bottom ******************************* 
        DO IPHASE=1,NPHASE
           K=2 ! Bottom
           S_VEL(:,:,K,IPHASE) = S_VEL(:,:,K,IPHASE) + REAL( 1 - P_BCS_ON(:,:,K-1,IPHASE) ) * INFINY*VEL_BCS(:,:,K,IPHASE) ! Dirichlet bc
           H_VEL(:,:,K,IPHASE,IPHASE) = REAL( 1 - P_BCS_ON(:,:,K-1,IPHASE) ) *H_VEL(:,:,K,IPHASE,IPHASE) + INFINY ! Dirichlet bc
           K=NZ ! Top
           S_VEL(:,:,K,IPHASE) = S_VEL(:,:,K,IPHASE) + REAL( 1 - P_BCS_ON(:,:,K,IPHASE) ) *INFINY*VEL_BCS(:,:,K,IPHASE) ! Dirichlet bc
           H_VEL(:,:,K,IPHASE,IPHASE) = H_VEL(:,:,K,IPHASE,IPHASE) + REAL( 1 - P_BCS_ON(:,:,K,IPHASE) ) *INFINY*VEL_BCS(:,:,K,IPHASE) ! Dirichlet bc
        END DO
        DO ITS_NONLIN_VEL=1,MAX(1,NITS_NONLIN_VEL) ! Non-linear iteration for limiting
! The advection velocity matrix term..
           DO K=2,NZ
           DO IUP=1,2
              CONSERV_VERT_ADV_VEL(IUP,:,:,K,:) = DEN_WATER_NOD(:,:,K,:)*0.5*(VEL(:,:,K-2+IUP,:)+VEL(:,:,K-1+IUP,:))/DX(3)
           END DO
           END DO
! apply limiting
           IF(NITS_NONLIN_VEL.NE.0) CALL AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV_VEL, VEL, VEL, NX+1,NY+1,NZ+1,NPHASE)
! Solve for velocity...
           CALL DIFFUSION_ADV_CAL(VEL,VEL2,KAPPA_VEL,CONSERV_VERT_ADV_VEL,CONSERV_VERT_ADV_VEL,H_VEL,S_VEL, &
                    NX+1,NY+1,NZ+1,NPHASE, .FALSE., VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER, 1, 0.0)
        END DO ! DO ITS_NONLIN_VEL=1,MAX(1,NITS_NONLIN_VEL)


! *******************************
! Calculate pressure...

! Calculate pressure matrix with matrix colouring method...        
        ALLOCATE(PCOLOR(NX,NY,NZ),VEC_COL_LONG(NX,NY,NZ,NPHASE),VEC_COL(NX,NY,NZ,NCOLOR)) 
        ALLOCATE(KAPPA_VOLFRA(2,3,NX+1,NY+1,NZ,NPHASE) )
        ALLOCATE(CONSERV_VERT_ADV_VOLFRA(2,NX+1,NY+1,NZ,NPHASE), CONSERV_VERT_ADV_VOLFRA_SIGN(2,NX+1,NY+1,NZ,NPHASE) )
        DO ICOLOR=1,3

           PCOLOR=0.0
           DO K=1+ICOLOR,NZ-1,3
              PCOLOR(:,:,K)=1.0
           END DO
        
           VEL2=0.0
           DO K=2,NZ
           DO J=2,NY
           DO I=2,NX
              MAT(:,:) = H_VEL(I,J,K,:,:)
              DO IPHASE=1,NPHASE
                 MAT(IPHASE,IPHASE) = MAT(IPHASE,IPHASE) + DEN_WATER_NOD(I,J,K,IPHASE)/DT
              END DO ! DO IPHASE=1,NPHASE
!              MAT = INVERSE(MAT)
              CALL INVERS_LEGS(MAT,NPHASE)
! VEL2= pivit_matrix^{-1} * C * pcolor
              DO IPHASE=1,NPHASE
                 DO JPHASE=1,NPHASE
                    VEL2(I,J,K,IPHASE)=VEL2(I,J,K,IPHASE) + MAT(IPHASE,JPHASE)*(PCOLOR(I,J,K) - PCOLOR(I,J,K-1))/DX(3)
                 END DO ! DO JPHASE=1,NPHASE
              END DO ! DO IPHASE=1,NPHASE
              
           END DO ! DO JPHASE=1,NPHASE
           END DO
           END DO
           
! Matrix vector multiplication involving the VEL2 that is VEC_COL=C^T * VEL2
           CONSERV_VERT_ADV_VOLFRA=0.0; CONSERV_VERT_ADV_VOLFRA_SIGN=0.0
           DO K=2,NZ
           DO IUP=1,2
!              CONSERV_VERT_ADV_VOLFRA(IUP,:,:,K,:) = DEN_WATER(:,:,K-1+IUP,:)*VEL2(:,:,K-1+IUP,:)/DX(3)
              CONSERV_VERT_ADV_VOLFRA(IUP,:,:,K,:) = DEN_WATER(:,:,K-1+IUP,:)*0.5*(1.0+SIGN(1.0,VEL(:,:,K-1+IUP,:)))*VEL2(:,:,K-1+IUP,:)/DX(3)  &
                                                   + DEN_WATER(:,:,K  +IUP,:)*0.5*(1.0-SIGN(1.0,VEL(:,:,K  +IUP,:)))*VEL2(:,:,K-1+IUP,:)/DX(3)
              CONSERV_VERT_ADV_VOLFRA_SIGN(IUP,:,:,K,:) = VEL(:,:,K-1+IUP,:) ! Use the sign to work out where info is comming from. 
           END DO
           END DO
! apply limiting
           IF(NITS_NONLIN_VOLFRA.NE.0) CALL AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV_VOLFRA, VEL, VOLFRA*DEN_WATER, NX+1,NY+1,NZ,NPHASE)
           H_VOLFRA=0.0; S_VOLFRA=0.0; KAPPA_VOLFRA=0.0 ! Set to zero as we are only doing C^T matrix vector multiplication.
! matrix vector
           CALL DIFFUSION_ADV_CAL(VOLFRA,VEC_COL_LONG,KAPPA_VOLFRA,CONSERV_VERT_ADV_VOLFRA,CONSERV_VERT_ADV_VOLFRA_SIGN, H_VOLFRA,S_VOLFRA, &
                        NX+1,NY+1,NZ,NPHASE, .TRUE., VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER, 1, 0.0) ! Solver not used
           DO K=2,NZ-1
           DO J=2,NY-1
           DO I=2,NX-1
              VEC_COL(I,J,K,ICOLOR) = SUM(VEC_COL_LONG(I,J,K,:) /DEN_WATER(I,J,K,:))  &
! Put compressible term into matrix
                                    + SUM(D_DEN_D_P(I,J,K,:)*VOLFRA(I,J,K,:)/DEN_WATER(I,J,K,:)) *  PCOLOR(I,J,K)
           END DO 
           END DO
           END DO

        END DO ! DO ICOLOR=1,3


! Caclaulate CT_VEL = C^T * VEL...
! Matrix vector multiplication involving the VEL that is VEC_COL=C^T * VEL
        DO K=2,NZ-1
        DO IUP=1,2
           CONSERV_VERT_ADV_VOLFRA(IUP,:,:,K,:) = DEN_WATER(:,:,K-2+IUP,:)*MAX(VEL(:,:,K-1+IUP,:),0.0)/DX(3)  &
                                                + DEN_WATER(:,:,K-1+IUP,:)*MIN(VEL(:,:,K-1+IUP,:),0.0)/DX(3) 
        END DO
        END DO
        ALLOCATE(CT_VEL(NX+1,NY+1,NZ,NPHASE), DP(NX+1,NY+1,NZ)) 
! limiting...
        IF(NITS_NONLIN_VOLFRA.NE.0) CALL AMEND_FOR_LIMITING_CONSERV_VERT_ADV(CONSERV_VERT_ADV_VOLFRA, VEL, VOLFRA*DEN_WATER, NX+1,NY+1,NZ,NPHASE)
! matrix vector...
        CALL DIFFUSION_ADV_CAL(VOLFRA,CT_VEL,KAPPA_VOLFRA,CONSERV_VERT_ADV_VOLFRA,CONSERV_VERT_ADV_VOLFRA, H_VOLFRA,S_VOLFRA, &
                     NX+1,NY+1,NZ,NPHASE, .TRUE., VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, 1, 0.0) ! Solver not used

        
! SOLVE FOR PRESSURE (use Thomas tri-digonal solver)...
        DO K=2,NZ-1
        DO J=2,NY-1
        DO I=2,NX-1
          DP(I,J,K) = SUM(CT_VEL(I,J,K,:)/DEN_WATER(I,J,K,:)) ! This is the rhs of the eqn that is being solved 
! - DIVID by density of water which is the normalization for global cty eqn
        END DO
        END DO
        END DO


        ALLOCATE(ABC(NX,NY,NZ,3))
        DO ICOLOR=1,3

           DO K=1+ICOLOR,NZ-1,3
           DO JCOLOR=1,3
              KK=K-1+JCOLOR
              ABC(:,:,K+JCOLOR-2, 4-JCOLOR)=VEC_COL(:,:,KK, ICOLOR) ! Extract tri-diagonal matrix.
           END DO
           END DO

        END DO ! DO ICOLOR=1,3

        DO J=2,NY-1
        DO I=2,NX-1
           DP(I,J,1)=0.0; DP(I,J,NZ)=0.0 ! Set halos of pressure to zero then solve for pressure...
           CALL tridag( ABC(I,J,2:NZ-1,1), ABC(I,J,2:NZ-1,2), ABC(I,J,2:NZ-1,3), DP(I,J,2:NZ-1),NZ-2)
        END DO
        END DO
        P = P + DP ! Update pressure


! *******************************
! Calculate velocity correction...
        VEL=0.0
        DO K=2,NZ
        DO J=2,NY
        DO I=2,NX

           MAT(:,:) = H_VEL(I,J,K,:,:)
           DO IPHASE=1,NPHASE
              MAT(IPHASE,IPHASE) = MAT(IPHASE,IPHASE) + DEN_WATER_NOD(I,J,K,IPHASE)/DT
           END DO

           CALL INVERS_LEGS(MAT,NPHASE) ! Invert matrix MAT

! Treat interphase velocity coupling implicitly in pressure...
           DO IPHASE=1,NPHASE
              DO JPHASE=1,NPHASE
                 VEL(I,J,K,IPHASE)=VEL(I,J,K,IPHASE) + MAT(iphase,jphase)*(DP(I,J,K) - DP(I,J,K-1))/DX(3)
              END DO ! DO JPHASE=1,NPHASE
           END DO ! DO IPHASE=1,NPHASE
        END DO
        END DO
        END DO

        CALL DEN_FUNCTION_WATER(DEN_WATER,D_DEN_D_P, TW,P, NX,NY,NZ,NPHASE) ! Use steam tables from FETCH

        RETURN
        END SUBROUTINE FLUIDS_CALC




        SUBROUTINE DEN_FUNCTION_WATER(DEN_WATER,D_DEN_D_P, TW,P, NX,NY,NZ,NPHASE)
! Use steam tables from FETCH to calculate DEN_WATER,D_DEN_D_P
! DEN_WATER=water density of all the phases
! D_DEN_D_P=derivative of density w.r.t. pressure.
         IMPLICIT NONE
         INTEGER, intent( in ) :: NX,NY,NZ,NPHASE
         REAL, intent( in ) :: DEN_WATER(NX+1,NY+1,NZ,NPHASE), TW(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( inout ) :: D_DEN_D_P(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
! Local variables...
         REAL, DIMENSION( : , : , :, : ), allocatable :: VEL_CV
         INTEGER :: I,J,K,IPHASE

        RETURN
        END SUBROUTINE DEN_FUNCTION_WATER



        SUBROUTINE CALC_VOLFRA_MASS_EXCHANGE(H_VOLFRA, DEN_WATER, VOLFRA, VEL, P, NX,NY,NZ,NPHASE) 
! Mass exchange term between the phases in H_VOLFRA
         IMPLICIT NONE
         INTEGER, intent( in ) :: NX,NY,NZ,NPHASE
         REAL, intent( inout ) :: H_VOLFRA(NX+1,NY+1,NZ,NPHASE,NPHASE) ! DEFINES MASS EXCHANGE
         REAL, intent( in ) :: VOLFRA(NX+1,NY+1,NZ,NPHASE)
         REAL, intent( in ) :: DEN_WATER(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
         REAL, intent( in ) :: VEL(NX+1,NY+1,NZ+1,NPHASE)
! Local variables...
         REAL, DIMENSION( : , : , :, : ), allocatable :: VEL_CV
         INTEGER :: I,J,K,IPHASE

         ALLOCATE(VEL_CV(NX+1, NY+1, NZ, NPHASE))
         DO K=1,NZ
            VEL_CV(:,:,K,:) = 0.5*( VEL(:,:,K,:) + VEL(:,:,K+1,:) ) 
         END DO
        RETURN
        END SUBROUTINE CALC_VOLFRA_MASS_EXCHANGE




        SUBROUTINE CALC_VEL_ABSORPTION(H_VEL, DEN_WATER_NOD, VOLFRA_NOD, VEL, P, NX,NY,NZ,NPHASE) 
! Calculate friction with pipe and between phases in H_VEL.
         IMPLICIT NONE
         INTEGER, intent( in ) :: NX,NY,NZ,NPHASE
         REAL, intent( inout ) :: H_VEL(NX+1,NY+1,NZ+1,NPHASE,NPHASE)
         REAL, intent( in ) :: VOLFRA_NOD(NX+1,NY+1,NZ+1,NPHASE)
         REAL, intent( in ) :: DEN_WATER_NOD(NX+1,NY+1,NZ,NPHASE),P(NX+1,NY+1,NZ)
         REAL, intent( inout ) :: VEL(NX+1,NY+1,NZ+1,NPHASE)
        RETURN
        END SUBROUTINE CALC_VEL_ABSORPTION




      subroutine tridag(a,b,c,d,nn)
      implicit none
      integer :: nn
! Solution of tridiagonal systems of equations
! The Thomas Algorithm is a special form of Gauss elimination that can be used to solve tridiago-
! nal systems of equations. 
! In CFD methods this algorithm is usually coded directly into the solution procedure, unless
! machine optimized subroutines are employed on a specific computer. A sample FORTRAN pro-
! gram to implement this algorithm is given here as:
!     solves a tridiagonal system using the Thomas Algorithm
!     there are nn equations, in the tridiagonal form:
!     a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
!     here, a(1) and c(nn) are assumed 0, and ignored
!     x is returned in d, b is altered
!     code set up to run on WATFOR-77
!     w.h. mason, April 10, 1992
      real a(nn),b(nn),c(nn),d(nn)
! local variables...
      integer i, k, km1
       real xm
!      if(nn .eq. 1)         then ! removed nn=1 case for speed.
!                            d(1)=d(1)/b(1)
!                            return
!                            end if
      do k = 2,nn
         km1     = k - 1
!      if(b(k-1) .eq. 0.0)   then ! remove checking for speed.
!                            write(6,100) km1
!                            stop
!                            end if
         xm      = a(k)/b(km1)
         b(k)    = b(k) - xm*c(km1)
         d(k)    = d(k) - xm*d(km1)
      end do ! do k = 2,nn

      d(nn)   = d(nn)/b(nn)
      k       = nn
      do i = 2,nn
         k       = nn + 1 - i
         d(k)    = (d(k) - c(k)*d(k+1))/b(k)
      end do ! do i = 2,nn
      return
      end subroutine tridag




        SUBROUTINE AXIS_SYM_DIFF_FUEL_RODS(TF,TFOLD, TW,TWOLD, SIGMA_ROD,  DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NR,NPHASE, &
                                           S_FUEL_RODS, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER)
! Axis-symmetric diffusion of temperature for rods...
! TF,TFOLD are the lattest and previos time level values of rod temperatures. 
! TW,TWOLD are the lattest and previos time level values of water temperatures. 
! DEN_CP_FUEL_RODS =CP*ROW of the fuel rods (could be control rods and includes cladding materials). 
! ROD_RADIUS_NODES are the nodes of the cells (side of the cells) for the model of the fuel pin in the radial direction.
! SIGMA_ROD contains the b.c's that we use to apply the fluid temp TW b.c's from to the rod temp. 
! If one need a non-uniform mesh spacing in radial direction then put this adjustment into the diffusion coefficient 
! in the radial direction. 
!
! x-coord will be along x-direction. 
! y-coord will be along radial direction.  x-direction. 
! z-coord will be along z-direction. 
! g-coord will be along radial direction. 
! NR is the no of cells in the radial direction of the fuel rod
! NX,NY,NZ no of cells in the x,y,z directions and the no of fuel rods = NX*NY (THIS INCLUDES HALOS)
! SOLVER OPTIONS:
! The smaller TF_RELAX the better the convergence =0,1,2,3 
! TF_RELAX = 0 : tridiagonal Gauss-Seidel.  
! TF_RELAX = 1 : pt Gauss-Seidel.  
! TF_RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems. 
! TF_RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) 
!   (only use for small nos of energy groups<10 as we invert a NGxNG matrix). - great for interphase coupling. 
! TF_MAX_ITS = max no of iterations.
! TF_ERROR_TOLER = error tolerence for convergence mas=x difference between temps from 2 consecutive iterations. 
! Similarly for other solver options. 
         IMPLICIT NONE
         INTEGER, intent( in ) :: NX,NY,NZ,NR,NPHASE, TF_MAX_ITS
         REAL, intent( in ) :: DX(3),ROD_RADIUS_NODES(NR+1),DT
         REAL, intent( inout ) :: TF(NX,NY,NZ,NR),TFOLD(NX,NY,NZ,NR), SIGMA_ROD(NX,NY,NZ,NPHASE) 
! Its NPHASE in SIGMA_ROD as we have to distribute cooling from water/gas phases.
         REAL, intent( in ) :: TW(NX+1,NY+1,NZ,1),TWOLD(NX+1,NY+1,NZ,1)
         REAL, intent( in ) :: DEN_CP_FUEL_RODS(NX,NY,NZ,NR), DIFF_FUEL_RODS(NX,NY,NZ,NR), S_FUEL_RODS(NX,NY,NZ,NR)
         REAL, intent( in ) :: TF_ERROR_TOLER
         INTEGER, intent( in ) :: TF_RELAX
! Local variables...
         REAL, PARAMETER :: TOLER = 1.e-14 ! Used to ensure we dont divid by zero. 
         INTEGER :: NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP
         REAL, DIMENSION(:,:,:,:), allocatable :: TF_TEMP, TFOLD_TEMP, DEN_CP_FUEL_RODS_TEMP, S_TEMP, DIFF_FUEL_RODS_TEMP
         REAL, DIMENSION(:,:,:,:), allocatable :: TF2_TEMP
         REAL, DIMENSION(:,:,:,:), allocatable :: H_DIAG_TEMP
         REAL, DIMENSION(:), allocatable ::  MID_RADIUS
         REAL, DIMENSION(:,:,:,:,:,:), allocatable :: KAPPA_FUEL
         REAL, DIMENSION(:,:, :,:, :), allocatable :: CONSERV_VERT_ADV(:,:, :,:, :)
         INTEGER :: I,J,K,G,IPHASE

! tranform to the following coords FOR SOLUTION:
! x-coord will be along radial direction. 
! y-coord will be along y-direction.  
! z-coord will be along x-direction. 
! g-coord will be along radial direction. 
        NX_TEMP=NR ! NO of radial cells
        NY_TEMP=NY
        NZ_TEMP=NZ
        NG_TEMP=NX
! 
! Interchange X and Z coords
        ALLOCATE(TF_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP), TFOLD_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP))
        ALLOCATE(DEN_CP_FUEL_RODS_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP), DIFF_FUEL_RODS_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP)) 
        DO I=1,NX
        DO G=1,NR
           TF_TEMP(G,:,:,I)=TF(I,:,:,G)
           TFOLD_TEMP(G,:,:,I)=TFOLD(I,:,:,G)
           DEN_CP_FUEL_RODS_TEMP(G,:,:,I) = DEN_CP_FUEL_RODS(I,:,:,G)
           DIFF_FUEL_RODS_TEMP(G,:,:,I) = DIFF_FUEL_RODS(I,:,:,G)
        END DO
        END DO   

! The radius between CVs 'nodes':
        ALLOCATE(MID_RADIUS(NR))
        DO I=1,NR
           MID_RADIUS(I)=0.5*(ROD_RADIUS_NODES(I)+ROD_RADIUS_NODES(I+1))
        END DO   
!
! calculate H_S_KAPPA for the fuel rods...
! We have to integrate over volume with this... 
        ALLOCATE( H_DIAG_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP),  S_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP) )
        H_DIAG_TEMP=0.0
        S_TEMP=0.0
        DO I=1,NX
        DO G=1,NR
           H_DIAG_TEMP(G,:,:,I)= MID_RADIUS(G) * DEN_CP_FUEL_RODS_TEMP(G,:,:,I) / DT 
           S_TEMP(G,:,:,I)= MID_RADIUS(G)*( S_FUEL_RODS(I,:,:,G) + DEN_CP_FUEL_RODS_TEMP(G,:,:,I)*TFOLD_TEMP(G,:,:,I)/DT  )
        END DO
        END DO
!
! Diffusion coefficient amended by radius...
! Set zero flux at the top of the fuel and bottom and at max radius of fuel rod
        ALLOCATE(KAPPA_FUEL(2,3,NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP))
        KAPPA_FUEL=0.0
        DO G=1,NG_TEMP ! No of energy groups
        DO K=2,NZ_TEMP-1
        DO J=2,NY_TEMP-1
        DO I=2,NX_TEMP-1
           KAPPA_FUEL(1,1,I,J,K,G)=ROD_RADIUS_NODES(I)  *0.5*( DIFF_FUEL_RODS_TEMP(I-1,J,K,G)+DIFF_FUEL_RODS_TEMP(I,J,K,G) ) &
                                  /MAX( ((ROD_RADIUS_NODES(I+1)-ROD_RADIUS_NODES(I))*(MID_RADIUS(I)-MID_RADIUS(I-1))), TOLER)
           KAPPA_FUEL(2,1,I,J,K,G)=ROD_RADIUS_NODES(I+1)*0.5*( DIFF_FUEL_RODS_TEMP(I+1,J,K,G)+DIFF_FUEL_RODS_TEMP(I,J,K,G) ) &
                                  /MAX( ((ROD_RADIUS_NODES(I+1)-ROD_RADIUS_NODES(I))*(MID_RADIUS(I+1)-MID_RADIUS(I))), TOLER)

           KAPPA_FUEL(1,3,I,J,K,G)=MID_RADIUS(I) *0.5*( DIFF_FUEL_RODS_TEMP(I,J,K-1,G)+DIFF_FUEL_RODS_TEMP(I,J,K,G) ) /DX(3)**2 ! Z is still Z thus using DX(3)
           KAPPA_FUEL(2,3,I,J,K,G)=MID_RADIUS(I) *0.5*( DIFF_FUEL_RODS_TEMP(I,J,K+1,G)+DIFF_FUEL_RODS_TEMP(I,J,K,G) ) /DX(3)**2
        END DO
        END DO
        END DO
        END DO
! Add zero flux bcs to KAPPA_FUEL (top, bottom, side of fuel/control rods):
        K=2 ! Bottom bc...
        KAPPA_FUEL(1,3,:,:,K,:)=0.0
        K=NZ_TEMP-1 ! Top bc...
        KAPPA_FUEL(2,3,:,:,K,:)=0.0
        I=NR-1 ! Side bc...
        KAPPA_FUEL(2,1,I,:,:,:)=0.0
!
! Set the temp b.c. of the fuel rods through SIGMA_ROD and using TW ******
        G=NR-1 ! The last row is where we set the b.c. next to the max radius of rod. 
        DO IPHASE=1,NPHASE
        DO I=2,NX-1
           S_TEMP(G,:,:,I)=  S_TEMP(G,:,:,I) + SIGMA_ROD(I,:,:,IPHASE) * TW(I,:,:,IPHASE)
           H_DIAG_TEMP(G,:,:,I)=  H_DIAG_TEMP(G,:,:,I) + SIGMA_ROD(I,:,:,IPHASE) 
        END DO 
        END DO 

         
        ALLOCATE(CONSERV_VERT_ADV(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP,2)); CONSERV_VERT_ADV=0.0
! CALCULATE TEMP OF CONTROL RODS
        ALLOCATE(TF2_TEMP(NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP))
        CALL DIFFUSION_ADV_CAL_DIAG_H(TF_TEMP,TF2_TEMP,KAPPA_FUEL,CONSERV_VERT_ADV,CONSERV_VERT_ADV,H_DIAG_TEMP,S_TEMP, &
                    NX_TEMP,NY_TEMP,NZ_TEMP,NG_TEMP, .FALSE., TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, 1, 0.0)

! Swap solution back - interchange X and Z coords
        DO I=1,NX
        DO G=1,NR
! ORIGINALY EXCHANGE:  TF_TEMP(G,:,:,I)=TF(I,:,:,G)
           TF(I,:,:,G)=TF_TEMP(G,:,:,I)
        END DO
        END DO     

        RETURN
        END SUBROUTINE AXIS_SYM_DIFF_FUEL_RODS




        SUBROUTINE NEUTRONICS_DIFFUSION(PSI,PSIOLD,S_NEUTRONS,SPEED,NEUTRON_DIFF,SIGMA_MATRIX,DT,DX, NX,NY,NZ,NG, &
                                        NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER,NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G, &
                                        NDELAY,BETA_EFF,BETA,LAMBDA,XI_D,XI_P,SIGMA_F_MATRIX,SIGMA_F,C) 
! This sub solves the neutron multi-energy group diffusion equation for the scalar flux PSI.
! PSIOLD is the old time level scalar flux. 
! C is delayed neutron precursor concentration. 
! NX is the no of cells in the x-direction including halo cells. 
! NY is the no of cells in the y-direction including halo cells. 
! NZ is the no of cells in the z-direction including halo cells. 
! NG is the no of energy groups.  
! PSI contains the best initial guess for the solver or the scalar flux for the matrix vector. 
! KAPPA(I,J,K,G,1or2,1or2or3) contains the diffusion coefficient between the control volume boundaries 
! for cell I,J,K and group G; the 1o2 indicates if its the left or right or bottom or top boundary;
! the 1or2or3 is the coordinate X,Y or Z. 
! dx(1:3) are the cell dimesions in all 3 directions
! dt is the time step size.
! H is the scatter removal operator. S_NEUTRONS is the source of neutrons
! S  is the source for the solver.
! SOLVER OPTIONS: 
! The smaller NEU_RELAX the better the convergence =0,1,2,3 
! NEU_RELAX = 0 : tridiagonal Gauss-Seidel.  
! NEU_RELAX = 1 : pt Gauss-Seidel.  
! NEU_RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems. 
! NEU_RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) (only use for small nos of energy groups<10 as we invert a NGxNG matrix).
! NEU_MAX_ITS = max no of iterations.
! NEU_ERROR_TOLER= max difference between the scalar flux between 2 consecutive iterations before assuming convergence. 
! NEU_MAX_ITS_IN_G = max no of iterations for the inner within group solve.
! NEU_ERROR_TOLER_IN_G = error tolerence for the inner within group solve.
! SOURCES AND PROBLEM SET UP IN S_NEUTRONS...
       IMPLICIT NONE
       INTEGER, intent( in ) :: NX,NY,NZ,NG, NEU_MAX_ITS, NEU_MAX_ITS_IN_G
       REAL, intent( inout ) :: PSI(NX,NY,NZ,NG)
       REAL, intent( in ) :: PSIOLD(NX,NY,NZ,NG), S_NEUTRONS(NX,NY,NZ,NG)
       REAL, intent( in ) :: NEUTRON_DIFF(NX,NY,NZ,NG)
       REAL, intent( in ) :: DX(3), NEU_ERROR_TOLER,NEU_ERROR_TOLER_IN_G, SPEED(NG), DT
       INTEGER, intent( in ) :: NEU_RELAX
       INTEGER, intent( in ) :: NDELAY
       REAL, intent( in ) :: BETA_EFF,BETA(NDELAY),LAMBDA(NDELAY),XI_D(NX,NY,NZ,NG),XI_P(NX,NY,NZ,NG)
       REAL, intent( in ) :: SIGMA_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F(NX,NY,NZ,NG),C(NX,NY,NZ,NDELAY)
! Local variables...
       REAL, DIMENSION( :,:, :,: ), allocatable :: PSI2, S
       REAL, DIMENSION( :,:, :,:, : ), allocatable :: H
       REAL, DIMENSION( :,:, :,:, :,: ), allocatable :: KAPPA
       REAL, DIMENSION( :,:, :,:, : ), allocatable :: CONSERV_VERT_ADV
       LOGICAL :: MATRIX_VEC

       ALLOCATE( PSI2(NX,NY,NZ,NG), KAPPA(2,3,NX,NY,NZ,NG), H(NX,NY,NZ,NG,NG),S(NX,NY,NZ,NG) ) 

       CALL CALC_H_S_KAPPA_FROM_SIGMA(KAPPA,H,S,PSI, DX,DT, NX,NY,NZ,NG, SPEED, NEUTRON_DIFF, SIGMA_MATRIX, S_NEUTRONS, &
                                      NDELAY,BETA_EFF,LAMBDA,BETA,XI_P,XI_D,SIGMA_F_MATRIX,SIGMA_F,C)

       MATRIX_VEC=.FALSE.
       ALLOCATE(CONSERV_VERT_ADV(2,NX,NY,NZ,NG)); CONSERV_VERT_ADV=0.0
       CALL DIFFUSION_ADV_CAL(PSI,PSI2,KAPPA,CONSERV_VERT_ADV,CONSERV_VERT_ADV, H,S, &
                     NX,NY,NZ,NG, MATRIX_VEC, NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER, NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G)
       RETURN
       END SUBROUTINE NEUTRONICS_DIFFUSION



       
       SUBROUTINE CALC_H_S_KAPPA_FROM_SIGMA(KAPPA,H,S,PSI, DX,DT, NX,NY,NZ,NG, SPEED, NEUTRON_DIFF, SIGMA_MATRIX, S_NEUTRONS, &
                                            NDELAY,BETA_EFF,LAMBDA,BETA,XI_P,XI_D,SIGMA_F_MATRIX,SIGMA_F,C)
       IMPLICIT NONE
! This sub returns KAPPA,H,S for use in the solver and calculated from other variables passed down. 
! Assume a cell with NEUTRON_DIFF(I,J,K,G) <= 0.0 is outside the calculations domain
! =0 assume a corresponding vacume bc <0 assume a symmetry b.c.
! Set PSI=0.0 in cells where we are using vacume b.c. with as they are used for extrapolating 0 b.c's. 
       INTEGER, intent( in ) :: NX,NY,NZ,NG,NDELAY
       REAL, intent( in ) :: DX(3),DT
       REAL, intent( in ) :: SPEED(NG)
       REAL, intent( inout ) :: KAPPA(2,3,NX,NY,NZ,NG), H(NX,NY,NZ,NG,NG), S(NX,NY,NZ,NG), PSI(NX,NY,NZ,NG)
       REAL, intent( in ) :: NEUTRON_DIFF(NX,NY,NZ,NG), S_NEUTRONS(NX,NY,NZ,NG)
       REAL, intent( in ) :: BETA_EFF,LAMBDA(NDELAY),BETA(NDELAY),XI_P(NX,NY,NZ,NG),XI_D(NX,NY,NZ,NG)
       REAL, intent( in ) :: SIGMA_F_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F(NX,NY,NZ,NG),C(NX,NY,NZ,NDELAY)
       REAL, intent( in ) :: SIGMA_MATRIX(NX,NY,NZ,NG,NG)
! local variables...
       INTEGER :: I,J,K,G,GD,II,III,IDIM,D
       REAL :: D_EXTRAP, DIFF_EFFECTIVE_DIV_DX2, R_DIAG_DIFF,R_OFF_DIAG_DIFF, R_NOT_ZERO_DIFF_ON_BOUNDARY, R_SWITCH_BC_ON

! ****************Calculate source S: 
       DO GD=1,NG
       DO G=1,NG
         S(:,:,:,G)=S_NEUTRONS(:,:,:,G) + PSI(:,:,:,G)/(DT*SPEED(G)) 
         DO D=1,NDELAY
            S(:,:,:,G)=S(:,:,:,G) + LAMBDA(D)*XI_D(:,:,:,G)*SIGMA_F(:,:,:,G) * C(:,:,:,D)
         END DO
       END DO
       END DO

! ***************************************************************************************
! ***********diffusion boundary conditions -symmetry and extrapolation to zero flux *****
! ***********AND volume diffusion calculation *******************************************
! the extrapolation distance d_extrap = 0.7104/SIGMA_TRANSPORT
! SIGMA_TRANSPORT is the transport cross section. 
! The diffusion coefficient DIFF = 1/(3*SIGMA_TRANSPORT)
! D_EXTRAP = 0.7104 *3.0 * NEUTRON_DIFF. 
! the effective diffusion coefficient is D_EXTRAP on boundary. 

! Assume a cell with NEUTRON_DIFF(I,J,K,G) <= 0.0 is outside the calculations domain
! =0 assume a corresponding vacume bc <0 assume a symmetry b.c.
       KAPPA = 0.0
       DO G=1,NG ! No of energy groups
       DO K=2,NZ-1
       DO I=2,NX-1
       DO J=2,NY-1
! Use SIGN(A, B): The kind of the return value is that of A and B. If B\ge 0 then the result is ABS(A), else it is -ABS(A). 
! Vacume(=value) OR symmetry(=0.0) bc:
          R_DIAG_DIFF = NEUTRON_DIFF(I,J,K,G)
          D_EXTRAP = 0.7104*3.0* R_DIAG_DIFF ! extrapolation distance.
          DIFF_EFFECTIVE_DIV_DX2= 1.0/MAX(1.E-10, D_EXTRAP) 
! x-direction:
          DO II=1,2
             III=(II-1)*2 -1

! x-direction:
             IDIM=1
             R_OFF_DIAG_DIFF = NEUTRON_DIFF(I+III,J,K,G)
! repeated...
             R_NOT_ZERO_DIFF_ON_BOUNDARY = 1.0  -  0.5*(-SIGN(1.0,R_OFF_DIAG_DIFF)+1.0)
             R_SWITCH_BC_ON = 0.5*(-SIGN(1.0,-R_OFF_DIAG_DIFF)+1.0) 
             KAPPA(II,IDIM,I,J,K,G)= (  DIFF_EFFECTIVE_DIV_DX2  * R_SWITCH_BC_ON &
                                +(0.5*( R_DIAG_DIFF + R_OFF_DIAG_DIFF ) /DX(IDIM)**2) *(1.0 - R_SWITCH_BC_ON) &
                              ) * R_NOT_ZERO_DIFF_ON_BOUNDARY

! y-direction:
             IDIM=2
             R_OFF_DIAG_DIFF = NEUTRON_DIFF(I,J+III,K,G)
! repeated...
             R_NOT_ZERO_DIFF_ON_BOUNDARY = 1.0  -  0.5*(-SIGN(1.0,R_OFF_DIAG_DIFF)+1.0)
             R_SWITCH_BC_ON = 0.5*(-SIGN(1.0,-R_OFF_DIAG_DIFF)+1.0) 
             KAPPA(II,IDIM,I,J,K,G)= (  DIFF_EFFECTIVE_DIV_DX2  * R_SWITCH_BC_ON &
                                +(0.5*( R_DIAG_DIFF + R_OFF_DIAG_DIFF ) /DX(IDIM)**2) *(1.0 - R_SWITCH_BC_ON) &
                              ) * R_NOT_ZERO_DIFF_ON_BOUNDARY

! z-direction:
             IDIM=3
             R_OFF_DIAG_DIFF = NEUTRON_DIFF(I,J,K+III,G)
! repeated...
             R_NOT_ZERO_DIFF_ON_BOUNDARY = 1.0  -  0.5*(-SIGN(1.0,R_OFF_DIAG_DIFF)+1.0)
             R_SWITCH_BC_ON = 0.5*(-SIGN(1.0,-R_OFF_DIAG_DIFF)+1.0) 
             KAPPA(II,IDIM,I,J,K,G)= (  DIFF_EFFECTIVE_DIV_DX2  * R_SWITCH_BC_ON &
                                +(0.5*( R_DIAG_DIFF + R_OFF_DIAG_DIFF ) /DX(IDIM)**2) *(1.0 - R_SWITCH_BC_ON) &
                              ) * R_NOT_ZERO_DIFF_ON_BOUNDARY
          END DO ! DO II=1,2
       END DO
       END DO
       END DO
       END DO   
! KAPPA coule be -ve if we have lots of symmetry cells thus make sure its >=0.0
       KAPPA=MAX(KAPPA,0.0) 
! Set PSI=0.0 in cells where we are using vacume b.c. with as they are used for extrapolating 0 b.c's. 
! Set PSI=0.0 if -ve NEUTRON_DIFF; Set PSI=0.0 if NEUTRON_DIFF = 0.0 ...
       PSI(:,:,:,:) = PSI(:,:,:,:) * 0.5*(-SIGN(1.0,-NEUTRON_DIFF(:,:,:,:))+1.0) 

! ***************************H calc: 
       DO GD=1,NG
       DO G=1,NG
          H(:,:,:,G,GD)= SIGMA_MATRIX(:,:,:,G,GD) + SIGMA_F_MATRIX(:,:,:,G,GD)
       END DO 
       END DO

       DO G=1,NG
         H(:,:,:,G,G)= H(:,:,:,G,G) + 1.0/(DT*SPEED(G))  
       END DO

       RETURN
       END SUBROUTINE CALC_H_S_KAPPA_FROM_SIGMA





       SUBROUTINE DIFFUSION_ADV_CAL(PSI,PSI2,KAPPA,CONSERV_VERT_ADV,CONSERV_VERT_ADV_SIGN,H,S, & 
                                    NX,NY,NZ,NG, MATRIX_VEC, RELAX,MAX_ITS,ERROR_TOLER, MAX_ITS_IN_G,ERROR_TOLER_IN_G) 
! NX is the no of cells in the x-direction including halo cells. 
! NY is the no of cells in the y-direction including halo cells. 
! NZ is the no of cells in the z-direction including halo cells. 
! NG is the no of energy groups.  
! PSI2=temp work space for solver or result of matrix vector multiplication if MATRIX_VEC=.true.
! PSI contains the best initial guess for the solver or the scalar flux for the matrix vector. 
! KAPPA(I,J,K,G,1or2,1or2or3) contains the diffusion coefficient between the control volume boundaries 
! for cell I,J,K and group G; the 1o2 indicates if its the left or right or bottom or top boundary;
! the 1or2or3 is the coordinate X,Y or Z. 
! Velocity advection is also placed into KAPPA.
! Conservative advection velocity is in CONSERV_VERT_ADV and the sign of the advection vel for 
! upwinding is in CONSERV_VERT_ADV_SIGN
! H is the scatter removal operator. S is the source of neutrons
! SOLVER OPTIONS: 
! RELAX=relax the matrix solver or not if true get safer convergence
! MAX_ITS=max no of its
! ERROR_TOLER= max difference between the scalar flux between 2 consecutive iterations before assuming convergence. 
! IF MAX_ITS_IN_G.NE.0 THEN USE A TRI-DIAGONAL SWEEP SOLVER.
! MAX_ITS_IN_G is the maximum no of inner iterations for a given group and 
! ERROR_TOLER_IN_G is the error tolerence for this inner iteration. 
! The smaller RELAX the better the convergence =0,1,2,3 
! RELAX = 0 : tridiagonal Gauss-Seidel.  
! RELAX = 1 : pt Gauss-Seidel.  
! RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems.  
! RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) (only use for small nos of energy groups<10 as we invert a NGxNG matrix). 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NX,NY,NZ,NG, MAX_ITS, MAX_ITS_IN_G
      REAL, intent( inout ) :: PSI(NX,NY,NZ,NG), PSI2(NX,NY,NZ,NG)
      REAL, intent( in ) :: KAPPA(2,3,NX,NY,NZ,NG)
      REAL, intent( in ) :: CONSERV_VERT_ADV(2,NX,NY,NZ+1,NG), CONSERV_VERT_ADV_SIGN(2,NX,NY,NZ+1,NG)
      REAL, intent( in ) :: H(NX,NY,NZ,NG,NG),S(NX,NY,NZ,NG), ERROR_TOLER, ERROR_TOLER_IN_G
      INTEGER, intent( in ) :: RELAX
      LOGICAL, intent( in ) :: MATRIX_VEC
! local variables...
      REAL :: ERROR
      INTEGER :: I,J,K,G,IUP, II,JJ,KK, ITS, ITS_G
      LOGICAL :: TRI_DIAG_SWEEPS, ENERGY_IMPLICIT
      REAL :: R_RELAX
      REAL, DIMENSION(:,:, :,:, :,:), allocatable :: MATRIX_NO_H
      REAL, DIMENSION(:,:, :,:), allocatable :: MATRIX_NO_H_DIAG, DIAG
      REAL, DIMENSION(:,:, :), allocatable :: PSI_IN_G
      REAL, DIMENSION(:,:, :), allocatable :: RHSX(:),RHSY(:),RHSZ(:)
      REAL, DIMENSION(:,:), allocatable :: MAT, MATINV
      REAL, DIMENSION(:), allocatable :: VECG, VECG2
      REAL, DIMENSION(:,:, :,:, :), allocatable :: MATINV_STORE_BL_DIA

! Map to the matrix:
       ALLOCATE( MATRIX_NO_H(2,3,NX,NY,NZ,NG) )
       ALLOCATE( MATRIX_NO_H_DIAG(NX,NY,NZ,NG) )


       IF(RELAX==2) THEN ! Over relax the diagonal for stability
          R_RELAX=1.0
       ELSE ! just the diagonal
          R_RELAX=0.0
       ENDIF

! form more efficient matrix...
       DO G=1,NG ! No of energy groups

       DO K=2,NZ-1
       DO J=2,NY-1
       DO I=2,NX-1
! Matrix vector...
          MATRIX_NO_H(1:2,1:3,I,J,K,G) = - KAPPA(1:2,1:3,I,J,K,G) 
          MATRIX_NO_H(1,3,I,J,K,G) = MATRIX_NO_H(1,3,I,J,K,G) +   0.5*(1.0+SIGN(1.0,CONSERV_VERT_ADV_SIGN(1,I,J,K,G))) * CONSERV_VERT_ADV(1,I,J,K,G)
          MATRIX_NO_H(2,3,I,J,K,G) = MATRIX_NO_H(2,3,I,J,K,G) +   0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(2,I,J,K,G))) * CONSERV_VERT_ADV(2,I,J,K,G)
          MATRIX_NO_H_DIAG(I,J,K,G) = SUM( KAPPA(1:2,1:3,I,J,K,G) ) + 0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(1,I,J,K,G))) * CONSERV_VERT_ADV(1,I,J,K,G) &
                                                                    + 0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(2,I,J,K,G))) * CONSERV_VERT_ADV(2,I,J,K,G)
!          MATRIX_NO_H(1,3,I,J,K,G) = MATRIX_NO_H(1,3,I,J,K,G) +   MAX(CONSERV_VERT_ADV(1,I,J,K,G),0.0)
!          MATRIX_NO_H(2,3,I,J,K,G) = MATRIX_NO_H(2,3,I,J,K,G) +   MIN(CONSERV_VERT_ADV(2,I,J,K,G),0.0)
!          MATRIX_NO_H_DIAG(I,J,K,G) = SUM( KAPPA(1:2,1:3,I,J,K,G) ) + MIN(CONSERV_VERT_ADV(1,I,J,K,G),0.0)+ MAX(CONSERV_VERT_ADV(2,I,J,K,G),0.0)
       END DO
       END DO
       END DO

       END DO    

! 
       IF(MATRIX_VEC) THEN

       PSI2=0.0

       DO G=1,NG ! No of energy groups

       DO K=2,NZ-1
       DO J=2,NY-1
       DO I=2,NX-1
! Matrix vector...
          PSI2(I,J,K,G)=PSI2(I,J,K,G)  + SUM(H(I,J,K,G,:)*PSI(I,J,K,:))  &
                      +  MATRIX_NO_H_DIAG(I,J,K,G)*PSI(I,J,K,G) & ! DIAGONAL DIFFUSION   
                      +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)  & ! OFF DIAGONAL -X-DIRECTION      
                      +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)  & ! OFF DIAGONAL -Y-DIRECTION    
                      +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)   ! OFF DIAGONAL -Z-DIRECTION  
       END DO
       END DO
       END DO

       END DO   
                        

       ELSE ! for solving the equations...

       TRI_DIAG_SWEEPS = (RELAX==0) ! default tri-diagonal solver
       IF(TRI_DIAG_SWEEPS) THEN
          ALLOCATE(RHSX(NX),RHSY(NY),RHSZ(NZ))
       ELSE ! Gauss-Seidel solver
          ALLOCATE(PSI_IN_G(NX,NY,NZ))
       ENDIF

       ENERGY_IMPLICIT = (RELAX==3) ! Forward-Backward block Gauss-Seidel (Energy implicit)
       IF(ENERGY_IMPLICIT) THEN
          ALLOCATE(MAT(NG,NG), MATINV(NG,NG), VECG(NG), VECG2(NG),  &
                   MATINV_STORE_BL_DIA(NX,NY,NZ,NG,NG))
          DO K = 2,NZ-1 
          DO J = 2,NY-1 
          DO I = 2,NX-1 

             MAT(:,:) = H(I,J,K,:,:)
             DO G=1,NG 
                MAT(G,G) = MAT(G,G) + MATRIX_NO_H_DIAG(I,J,K,G)
             END DO
             MATINV=MAT
!             MATINV=invers(MATINV) 
             CALL INVERS_LEGS(MATINV,NG)
             MATINV_STORE_BL_DIA(I,J,K,:,:) = MATINV(:,:)

          END DO
          END DO
          END DO
          
       ENDIF

! Calculate diagonal...

       ALLOCATE(DIAG(NX,NY,NZ,NG)) 
       DO G=1,NG ! No of energy groups

       DO K=2,NZ-1
       DO J=2,NY-1
       DO I=2,NX-1
          DIAG(I,J,K,G) = MATRIX_NO_H_DIAG(I,J,K,G) + R_RELAX*SUM(ABS(H(I,J,K,G,:))) + (1.-R_RELAX)*H(I,J,K,G,G)
       END DO
       END DO
       END DO

       END DO

! Reactor calculation...
       
       DO ITS=1,MAX_ITS

       PSI2=PSI
       
       IF(TRI_DIAG_SWEEPS) THEN ! *****************
          DO G=1,NG ! No of energy groups

          DO ITS_G=1,MAX_ITS_IN_G

! X-direction
          IF(NX.GT.3) THEN ! Dont bother using a tridiagonal solver if only 1 point
          ERROR=0.0
          DO KK=0,1
          DO K=2*(1-KK)+ (NZ-1)*KK,  (NZ-1)*(1-KK)  +2*KK, 1*(1-KK)  -1*KK
          DO JJ=0,1
          DO J=2*(1-JJ)+ (NY-1)*JJ,  (NY-1)*(1-JJ)  +2*JJ, 1*(1-JJ)  -1*JJ
          RHSX=0.0
          DO I=2,NX-1
             RHSX(I) = S(I,J,K,G) - SUM(PSI(I,J,K,:)*H(I,J,K,G,:)) + PSI(I,J,K,G)*H(I,J,K,G,G) & ! Only off diag of H 
                     -(  & ! DIAGONAL DIFFUSION   
!                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             I=2 ! Include the potential bcs:
             RHSX(I) = RHSX(I) -MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) 
             I=NX-1 ! Include the potential bcs:
             RHSX(I) = RHSX(I) -MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,1,2:NX-1,J,K,G), MATRIX_NO_H_DIAG(2:NX-1,J,K,G) + H(2:NX-1,J,K,G,G), MATRIX_NO_H(2,1,2:NX-1,J,K,G), RHSX(2:NX-1), NX-2)
             ERROR=MAX( MAXVAL(ABS(PSI(2:NX-1,J,K,G) - RHSX(2:NX-1)) ), ERROR)
             PSI(2:NX-1,J,K,G) = RHSX(2:NX-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NX.GT.3) THEN

! Y-direction
          IF(NY.GT.3) THEN
          ERROR=0.0
          DO KK=0,1
          DO K=2*(1-KK)+ (NZ-1)*KK,  (NZ-1)*(1-KK)  +2*KK, 1*(1-KK)  -1*KK
          DO II=0,1
          DO I=2*(1-II)+ (NX-1)*II,  (NX-1)*(1-II)  +2*II, 1*(1-II)  -1*II
          RHSY=0.0
          DO J=2,NY-1
             RHSY(J) = S(I,J,K,G) - SUM(PSI(I,J,K,:)*H(I,J,K,G,:)) + PSI(I,J,K,G)*H(I,J,K,G,G) & ! Only off diag of H
                     -(  & ! DIAGONAL DIFFUSION   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
!                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             J=2 ! Include the potential bcs:
             RHSY(J) = RHSY(J) -MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) 
             J=NY-1 ! Include the potential bcs:
             RHSY(J) = RHSY(J) -MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,2,I,2:NY-1,K,G), MATRIX_NO_H_DIAG(I,2:NY-1,K,G) + H(I,2:NY-1,K,G,G), MATRIX_NO_H(2,2,I,2:NY-1,K,G), RHSY(2:NY-1), NY-2)
             ERROR=MAX( MAXVAL(ABS(PSI(I,2:NY-1,K,G) - RHSY(2:NY-1)) ), ERROR)
             PSI(I,2:NY-1,K,G) = RHSY(2:NY-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NY.GT.3) THEN

! Z-direction
          IF(NZ.GT.3) THEN
          ERROR=0.0
          DO JJ=0,1
          DO J=2*(1-JJ)+ (NY-1)*JJ,  (NY-1)*(1-JJ)  +2*JJ, 1*(1-JJ)  -1*JJ
          DO II=0,1
          DO I=2*(1-II)+ (NX-1)*II,  (NX-1)*(1-II)  +2*II, 1*(1-II)  -1*II
          RHSZ=0.0
          DO K=2,NZ-1
             RHSZ(K) = S(I,J,K,G) - SUM(PSI(I,J,K,:)*H(I,J,K,G,:)) + PSI(I,J,K,G)*H(I,J,K,G,G) & ! Only off diag of H
                     -(  & ! DIAGONAL DIFFUSION   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
!                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             K=2 ! Include the potential bcs:
             RHSX(K) = RHSX(K) -MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) 
             K=NZ-1 ! Include the potential bcs:
             RHSX(K) = RHSX(K) -MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,3,I,J,2:NZ-1,G), MATRIX_NO_H_DIAG(I,J,2:NZ-1,G) + H(I,J,2:NZ-1,G,G), MATRIX_NO_H(2,3,I,J,2:NZ-1,G), RHSZ(2:NZ-1), NZ-2)
             ERROR=MAX( MAXVAL(ABS(PSI(I,J,2:NZ-1,G) - RHSZ(2:NZ-1)) ), ERROR)
             PSI(I,J,2:NZ-1,G) = RHSZ(2:NZ-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NZ.GT.3) THEN

          END DO ! DO ITS_G=1,MAX_ITS_IN_G

          END DO ! DO G=1,NG

       ELSE IF(ENERGY_IMPLICIT) THEN ! IF(TRI_DIAG_SWEEPS) THEN ******************
! Forward-Backward block Gauss-Seidel (Energy implicit)

          DO IUP=0,1 ! This is forward backward Gauss-Seidel...
          DO K = 2*(1-IUP)+ (NZ-1)*IUP,  (NZ-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP
          DO J = 2*(1-IUP)+ (NY-1)*IUP,  (NY-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP
          DO I = 2*(1-IUP)+ (NX-1)*IUP,  (NX-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP

             DO G=1,NG ! No of energy groups
             VECG(G)= S(I,J,K,G) & ! block diag implicit   - SUM(PSI(I,J,K,:)*H(I,J,K,G,:))  - MATRIX_NO_H_DIAG(I,J,K,G)*PSI(I,J,K,G)  & ! Subtract out block diagonal
                     -(  & !  diag implicit (diagonal diffusion)  MATRIX_NO_H_DIAG(I,J,K,G)*PSI(I,J,K,G)   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )
             END DO ! DO G=1,NG

             DO G=1,NG 
                PSI(I,J,K,G) = SUM(MATINV_STORE_BL_DIA(I,J,K,G,:)* VECG2(:)) 
             END DO
             
             PSI(I,J,K,:)=VECG2(:) 

          END DO
          END DO
          END DO
          END DO ! DO IUP=0,1

       ELSE ! IF(TRI_DIAG_SWEEPS) THEN ELSE IF(ENERGY_IMPLICIT) THEN ******************
! Forward-Backward pt Gauss-Seidel 
          DO G=1,NG ! No of energy groups

          DO ITS_G=1,MAX_ITS_IN_G

          PSI_IN_G(:,:,:)=PSI(:,:,:,G)

          DO IUP=0,1 ! This is forward backward Gauss-Seidel...
          DO K = 2*(1-IUP)+ (NZ-1)*IUP,  (NZ-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP
          DO J = 2*(1-IUP)+ (NY-1)*IUP,  (NY-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP
          DO I = 2*(1-IUP)+ (NX-1)*IUP,  (NX-1)*(1-IUP)  +2*IUP, 1*(1-IUP)  -1*IUP
            PSI(I,J,K,G)= (S(I,J,K,G) +DIAG(I,J,K,G)*PSI(I,J,K,G)- SUM(PSI(I,J,K,:)*H(I,J,K,G,:))  &
                     -( MATRIX_NO_H_DIAG(I,J,K,G)*PSI(I,J,K,G) & ! DIAGONAL DIFFUSION   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G) ) & ! OFF DIAGONAL -Z-DIRECTION 
                           )/DIAG(I,J,K,G)
          END DO
          END DO
          END DO
          END DO ! DO IUP=0,1

       ERROR=MAXVAL( ABS(PSI_IN_G(:,:,:)-PSI(:,:,:,G)) ) 
print *, 'error_g', error, 'its_g', its_g 
       IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE

          END DO ! DO ITS_G=1,MAX_ITS_IN_G

          END DO ! DO G=1,NG
       ENDIF ! ENDOF IF(TRI_DIAG_SWEEPS) THEN ELSE ********

       
       ERROR=MAXVAL(ABS(PSI2-PSI)) 
       IF(ERROR.LT.ERROR_TOLER) CYCLE

       END DO ! DO ITS=1,MAX_ITS

       ENDIF ! ENDOF IF(MATRIX_VEC.NE.0) THEN
       RETURN
       END SUBROUTINE DIFFUSION_ADV_CAL
       
       


       SUBROUTINE DIFFUSION_ADV_CAL_DIAG_H(PSI,PSI2,KAPPA,CONSERV_VERT_ADV,CONSERV_VERT_ADV_SIGN,H_DIAG,S, & 
                                    NX,NY,NZ,NG, MATRIX_VEC, RELAX,MAX_ITS,ERROR_TOLER, MAX_ITS_IN_G,ERROR_TOLER_IN_G) 
! NX is the no of cells in the x-direction including halo cells. 
! NY is the no of cells in the y-direction including halo cells. 
! NZ is the no of cells in the z-direction including halo cells. 
! NG is the no of energy groups.  
! PSI2=temp work space for solver or result of matrix vector multiplication if MATRIX_VEC=.true.
! PSI contains the best initial guess for the solver or the scalar flux for the matrix vector. 
! KAPPA(I,J,K,G,1or2,1or2or3) contains the diffusion coefficient between the control volume boundaries 
! for cell I,J,K and group G; the 1o2 indicates if its the left or right or bottom or top boundary;
! the 1or2or3 is the coordinate X,Y or Z. 
! Velocity advection is also placed into KAPPA.
! Conservative advection velocity is in CONSERV_VERT_ADV and the sign of the advection vel for 
! upwinding is in CONSERV_VERT_ADV_SIGN
! H is the scatter removal operator. S is the source of neutrons
! SOLVER OPTIONS: 
! RELAX=relax the matrix solver or not if true get safer convergence
! MAX_ITS=max no of its
! ERROR_TOLER= max difference between the scalar flux between 2 consecutive iterations before assuming convergence. 
! IF MAX_ITS_IN_G.NE.0 THEN USE A TRI-DIAGONAL SWEEP SOLVER.
! MAX_ITS_IN_G is the maximum no of inner iterations for a given group and 
! ERROR_TOLER_IN_G is the error tolerence for this inner iteration. 
! The smaller RELAX the better the convergence =0,1,2,3 
! RELAX = 0 : tridiagonal Gauss-Seidel.  
! RELAX = 1 : pt Gauss-Seidel.  
! RELAX = 2 : relaxed pt Gauss-Seidel - relaxed for better converence in difficult problems.  
! RELAX = 3 : Forward-Backward block Gauss-Seidel (Energy implicit) (only use for small nos of energy groups<10 as we invert a NGxNG matrix). 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NX,NY,NZ,NG, MAX_ITS, MAX_ITS_IN_G
      REAL, intent( inout ) :: PSI(NX,NY,NZ,NG), PSI2(NX,NY,NZ,NG)
      REAL, intent( in ) :: KAPPA(2,3,NX,NY,NZ,NG)
      REAL, intent( in ) :: CONSERV_VERT_ADV(2,NX,NY,NZ+1,NG), CONSERV_VERT_ADV_SIGN(2,NX,NY,NZ+1,NG)
      REAL, intent( in ) :: H_DIAG(NX,NY,NZ,NG),S(NX,NY,NZ,NG), ERROR_TOLER, ERROR_TOLER_IN_G
      INTEGER, intent( in ) :: RELAX
      LOGICAL, intent( in ) :: MATRIX_VEC
! local variables...
      REAL :: ERROR
      INTEGER :: I,J,K,G,IUP, II,JJ,KK, ITS, ITS_G
      REAL :: R_RELAX
      REAL, DIMENSION(:,:, :,:, :,:), allocatable :: MATRIX_NO_H
      REAL, DIMENSION(:,:, :,:), allocatable :: MATRIX_NO_H_DIAG, DIAG
      REAL, DIMENSION(:,:, :), allocatable :: RHSX(:),RHSY(:),RHSZ(:)
      REAL, DIMENSION(:,:), allocatable :: MAT, MATINV
      REAL, DIMENSION(:), allocatable :: VECG, VECG2
      REAL, DIMENSION(:,:, :,:, :), allocatable :: MATINV_STORE_BL_DIA

! Map to the matrix:
       ALLOCATE( MATRIX_NO_H(2,3,NX,NY,NZ,NG) )
       ALLOCATE( MATRIX_NO_H_DIAG(NX,NY,NZ,NG) )


       IF(RELAX==2) THEN ! Over relax the diagonal for stability
          R_RELAX=1.0
       ELSE ! just the diagonal
          R_RELAX=0.0
       ENDIF

! form more efficient matrix...
       DO G=1,NG ! No of energy groups

       DO K=2,NZ-1
       DO J=2,NY-1
       DO I=2,NX-1
! Matrix vector...
          MATRIX_NO_H(1:2,1:3,I,J,K,G) = - KAPPA(1:2,1:3,I,J,K,G) 
          MATRIX_NO_H(1,3,I,J,K,G) = MATRIX_NO_H(1,3,I,J,K,G) +   0.5*(1.0+SIGN(1.0,CONSERV_VERT_ADV_SIGN(1,I,J,K,G))) * CONSERV_VERT_ADV(1,I,J,K,G)
          MATRIX_NO_H(2,3,I,J,K,G) = MATRIX_NO_H(2,3,I,J,K,G) +   0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(2,I,J,K,G))) * CONSERV_VERT_ADV(2,I,J,K,G)
          MATRIX_NO_H_DIAG(I,J,K,G) = SUM( KAPPA(1:2,1:3,I,J,K,G) ) + 0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(1,I,J,K,G))) * CONSERV_VERT_ADV(1,I,J,K,G) &
                                                                    + 0.5*(1.0-SIGN(1.0,CONSERV_VERT_ADV_SIGN(2,I,J,K,G))) * CONSERV_VERT_ADV(2,I,J,K,G)
!          MATRIX_NO_H(1,3,I,J,K,G) = MATRIX_NO_H(1,3,I,J,K,G) +   MAX(CONSERV_VERT_ADV(1,I,J,K,G),0.0)
!          MATRIX_NO_H(2,3,I,J,K,G) = MATRIX_NO_H(2,3,I,J,K,G) +   MIN(CONSERV_VERT_ADV(2,I,J,K,G),0.0)
!          MATRIX_NO_H_DIAG(I,J,K,G) = SUM( KAPPA(1:2,1:3,I,J,K,G) ) + MIN(CONSERV_VERT_ADV(1,I,J,K,G),0.0)+ MAX(CONSERV_VERT_ADV(2,I,J,K,G),0.0)
       END DO
       END DO
       END DO

       END DO    

! 
       IF(MATRIX_VEC) THEN

       PSI2=0.0

       DO G=1,NG ! No of energy groups

       DO K=2,NZ-1
       DO J=2,NY-1
       DO I=2,NX-1
! Matrix vector...
          PSI2(I,J,K,G)=PSI2(I,J,K,G)  + H_DIAG(I,J,K,G)*PSI(I,J,K,G)  &
                      +  MATRIX_NO_H_DIAG(I,J,K,G)*PSI(I,J,K,G) & ! DIAGONAL DIFFUSION   
                      +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)  & ! OFF DIAGONAL -X-DIRECTION      
                      +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)  & ! OFF DIAGONAL -Y-DIRECTION    
                      +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)   ! OFF DIAGONAL -Z-DIRECTION  
       END DO
       END DO
       END DO

       END DO   
                        

       ELSE ! for solving the equations...

       ALLOCATE(RHSX(NX),RHSY(NY),RHSZ(NZ))

! Reactor calculation...
       
       DO ITS=1,MAX_ITS

       PSI2=PSI
       
          DO G=1,NG ! No of energy groups

          DO ITS_G=1,MAX_ITS_IN_G

! X-direction
          IF(NX.GT.3) THEN ! Dont bother using a tridiagonal solver if only 1 point
          ERROR=0.0
          DO KK=0,1
          DO K=2*(1-KK)+ (NZ-1)*KK,  (NZ-1)*(1-KK)  +2*KK, 1*(1-KK)  -1*KK
          DO JJ=0,1
          DO J=2*(1-JJ)+ (NY-1)*JJ,  (NY-1)*(1-JJ)  +2*JJ, 1*(1-JJ)  -1*JJ
          RHSX=0.0
          DO I=2,NX-1
             RHSX(I) = S(I,J,K,G)  &
                     -(  & ! DIAGONAL DIFFUSION   
!                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             I=2 ! Include the potential bcs:
             RHSX(I) = RHSX(I) -MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) 
             I=NX-1 ! Include the potential bcs:
             RHSX(I) = RHSX(I) -MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,1,2:NX-1,J,K,G), MATRIX_NO_H_DIAG(2:NX-1,J,K,G) + H_DIAG(2:NX-1,J,K,G), MATRIX_NO_H(2,1,2:NX-1,J,K,G), RHSX(2:NX-1), NX-2)
             ERROR=MAX( MAXVAL(ABS(PSI(2:NX-1,J,K,G) - RHSX(2:NX-1)) ), ERROR)
             PSI(2:NX-1,J,K,G) = RHSX(2:NX-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NX.GT.3) THEN

! Y-direction
          IF(NY.GT.3) THEN
          ERROR=0.0
          DO KK=0,1
          DO K=2*(1-KK)+ (NZ-1)*KK,  (NZ-1)*(1-KK)  +2*KK, 1*(1-KK)  -1*KK
          DO II=0,1
          DO I=2*(1-II)+ (NX-1)*II,  (NX-1)*(1-II)  +2*II, 1*(1-II)  -1*II
          RHSY=0.0
          DO J=2,NY-1
             RHSY(J) = S(I,J,K,G)   &
                     -(  & ! DIAGONAL DIFFUSION   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
!                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             J=2 ! Include the potential bcs:
             RHSY(J) = RHSY(J) -MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) 
             J=NY-1 ! Include the potential bcs:
             RHSY(J) = RHSY(J) -MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,2,I,2:NY-1,K,G), MATRIX_NO_H_DIAG(I,2:NY-1,K,G) + H_DIAG(I,2:NY-1,K,G), MATRIX_NO_H(2,2,I,2:NY-1,K,G), RHSY(2:NY-1), NY-2)
             ERROR=MAX( MAXVAL(ABS(PSI(I,2:NY-1,K,G) - RHSY(2:NY-1)) ), ERROR)
             PSI(I,2:NY-1,K,G) = RHSY(2:NY-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NY.GT.3) THEN

! Z-direction
          IF(NZ.GT.3) THEN
          ERROR=0.0
          DO JJ=0,1
          DO J=2*(1-JJ)+ (NY-1)*JJ,  (NY-1)*(1-JJ)  +2*JJ, 1*(1-JJ)  -1*JJ
          DO II=0,1
          DO I=2*(1-II)+ (NX-1)*II,  (NX-1)*(1-II)  +2*II, 1*(1-II)  -1*II
          RHSZ=0.0
          DO K=2,NZ-1
             RHSZ(K) = S(I,J,K,G)  &
                     -(  & ! DIAGONAL DIFFUSION   
                     +  MATRIX_NO_H(1,1,I,J,K,G)*PSI(I-1,J,K,G) + MATRIX_NO_H(2,1,I,J,K,G)*PSI(I+1,J,K,G)   & ! OFF DIAGONAL -X-DIRECTION      
                     +  MATRIX_NO_H(1,2,I,J,K,G)*PSI(I,J-1,K,G) + MATRIX_NO_H(2,2,I,J,K,G)*PSI(I,J+1,K,G)   & ! OFF DIAGONAL -Y-DIRECTION    
!                     +  MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) + MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G)  & ! OFF DIAGONAL -Z-DIRECTION 
                           )  
          END DO
             K=2 ! Include the potential bcs:
             RHSX(K) = RHSX(K) -MATRIX_NO_H(1,3,I,J,K,G)*PSI(I,J,K-1,G) 
             K=NZ-1 ! Include the potential bcs:
             RHSX(K) = RHSX(K) -MATRIX_NO_H(2,3,I,J,K,G)*PSI(I,J,K+1,G) 
! Tri-diagonal solver: 
             CALL tridag(MATRIX_NO_H(1,3,I,J,2:NZ-1,G), MATRIX_NO_H_DIAG(I,J,2:NZ-1,G) + H_DIAG(I,J,2:NZ-1,G), MATRIX_NO_H(2,3,I,J,2:NZ-1,G), RHSZ(2:NZ-1), NZ-2)
             ERROR=MAX( MAXVAL(ABS(PSI(I,J,2:NZ-1,G) - RHSZ(2:NZ-1)) ), ERROR)
             PSI(I,J,2:NZ-1,G) = RHSZ(2:NZ-1) 
          END DO
          END DO
          END DO
          END DO
          IF(ERROR.LT.ERROR_TOLER_IN_G) CYCLE
          ENDIF ! IF(NZ.GT.3) THEN

          END DO ! DO ITS_G=1,MAX_ITS_IN_G

          END DO ! DO G=1,NG


       
       ERROR=MAXVAL(ABS(PSI2-PSI)) 
       IF(ERROR.LT.ERROR_TOLER) CYCLE

       END DO ! DO ITS=1,MAX_ITS

       ENDIF ! ENDOF IF(MATRIX_VEC.NE.0) THEN
       RETURN
       END SUBROUTINE DIFFUSION_ADV_CAL_DIAG_H
       
       


      SUBROUTINE invers_LEGS(A,N)
! sub return the invers of the matrix A and returns in AINV
! A is overwritten

! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
!
      INTEGER I,J,jj,N,INDX(N)
      real A(N,N),B(N),X(N),AINV(N,N)

!      ewrite(3,*) MAT,N,B
!
!      ewrite(3,*) 'In LEGS'
!      do I=1,N
!      do J=1,N
!        A(I,J) = MAT((J-1)*N+I)
!       enddo
!      enddo
      CALL ELGS(A,N,INDX)

      do jj=1,n
      b=0.0
      b(jj)=1.0
!
      do  I = 1, N-1! Was loop 100
      do  J = I+1, N! Was loop 90
            B(INDX(J)) = B(INDX(J)) -A(INDX(J),I)*B(INDX(I))
      end do ! Was loop 90
      end do ! Was loop 100
!
      X(N) = B(INDX(N))/A(INDX(N),N)
      do  I = N-1, 1, -1! Was loop 200
        X(I) = B(INDX(I))
      do  J = I+1, N! Was loop 190
          X(I) = X(I)-A(INDX(I),J)*X(J)
      end do ! Was loop 190
          X(I) =  X(I)/A(INDX(I),I)
      end do ! Was loop 200

      ainv(:,jj)=b(:)

      end do ! do jj=1,n

      a=ainv
!
      RETURN
      END SUBROUTINE invers_LEGS


      SUBROUTINE LEGS(MAT,N,B,X,A,INDX)
!
! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
!
      INTEGER I,J,N,INDX(N)
      real A(N,N),B(N),X(N),MAT(N*N)

!      ewrite(3,*) MAT,N,B
!
!      ewrite(3,*) 'In LEGS'
      do I=1,N
      do J=1,N
        A(I,J) = MAT((J-1)*N+I)
       enddo
      enddo
      CALL ELGS(A,N,INDX)
!
      do  I = 1, N-1! Was loop 100
      do  J = I+1, N! Was loop 90
            B(INDX(J)) = B(INDX(J)) -A(INDX(J),I)*B(INDX(I))
      end do ! Was loop 90
      end do ! Was loop 100
!
      X(N) = B(INDX(N))/A(INDX(N),N)
      do  I = N-1, 1, -1! Was loop 200
        X(I) = B(INDX(I))
      do  J = I+1, N! Was loop 190
          X(I) = X(I)-A(INDX(I),J)*X(J)
      end do ! Was loop 190
          X(I) =  X(I)/A(INDX(I),I)
      end do ! Was loop 200
!
      RETURN
      END SUBROUTINE LEGS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ELGS(A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.
      INTEGER I,J,N,K,itmp,INDX(N)
      real A(N,N),C(N),C1,PI,PI1,pj
! Initialize the index
!      ewrite(3,*) 'IN ELGS'
      do     I = 1, N! Was loop 50
        INDX(I) = I
      end do ! Was loop 50
!
! Find the rescaling factors, one from each row
!
      do    I = 1, N! Was loop 100
          C1= 0.0
      do    J = 1, N! Was loop 90
            C1 = AMAX1(C1,ABS(A(I,J)))
!            ewrite(3,*) C1
      end do ! Was loop 90
          C(I) = C1
      end do ! Was loop 100
!
! Search the pivoting (largest) element from each column
!
      do    J = 1, N-1! Was loop 200
        PI1 = 0.0
      do    I = J, N! Was loop 150
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
      end do ! Was loop 150
!
! Interchange the rows via INDX(N) to record pivoting order
!
!        ewrite(3,*) indx
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
      do    I = J+1, N! Was loop 170
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      do    K = J+1, N! Was loop 160
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      end do ! Was loop 160
      end do ! Was loop 170
      end do ! Was loop 200
!
      RETURN
      END SUBROUTINE ELGS
!

      subroutine  write_tw_btm(NX, NY, NZ, NPHASE, TW)
         implicit none
         integer, intent(in) :: NX, NY, NZ, NPHASE
         real, intent(in)    :: TW(NX+1,NY+1,NZ,NPHASE) 
         integer i, j, k, iphase, io

         io = 143
         open(io, file='water-temp_btm.dat')
         k=1
         write(io,*) "# to plot this at a terminal type... "
         write(io,*) "# gnuplot "
         write(io,*) "# set pm3d map "
         write(io,*) "# splot './water-temp_btm.dat' us 1:2:4 with pm3d "
         write(io,*) "# x y z nPhase value"
         do i=1,NX
         do j=1,NY
         do iphase=1,NPHASE
             write(io,'(3i8,f16.4)') i, j, iphase, tw(i,j,k,iphase)        
         enddo
         enddo
         write(io,'(3i8,f16.4)')
         enddo
         close(io)
     end subroutine  write_tw_btm

      subroutine  write_tw_top(NX, NY, NZ, NPHASE, TW)
         implicit none
         integer, intent(in) :: NX, NY, NZ, NPHASE
         real, intent(in)    :: TW(NX+1,NY+1,NZ,NPHASE) 
         integer i, j, k, iphase,  io

         io = 143
         open(io, file='water-temp_top.dat')
         k=NZ
         write(io,*) "# to plot this at a terminal type... "
         write(io,*) "# gnuplot "
         write(io,*) "# set pm3d map "
         write(io,*) "# splot './water-temp_top.dat' us 1:2:4 with pm3d "
         write(io,*) "# "
         write(io,*) "# x y z nPhase value"
         do i=1,NX
         do j=1,NY
         do iphase=1,NPHASE
             write(io,'(3i8,f16.4)') i, j, iphase, tw(i,j,k,iphase)        
         enddo
         enddo
         write(io,'(3i8,f16.4)')
         enddo
         close(io)
     end subroutine  write_tw_top

     subroutine write_tf(NX, NY, NZ, NR, TF, ROD_RADIUS_NODES)
         implicit none
         integer, intent(in) :: NX, NY, NZ, NR
         real, intent(in)    :: TF(NX,NY,NZ,NR)
         real, intent(in)    :: ROD_RADIUS_NODES(NR+1)
         integer i, j, k, io, IR

         io = 143
         open(143, file='fuel-temp.dat')
         write(io,*) "# to plot this at a terminal type..."
         write(io,*) "# gnuplot "
         write(io,*) "# pl './fuel-temp.dat' w linesp pt 7 ps 1.2 "
         write(io,*) "#  "
         write(143,*) "# radial node , value for cell i=2,j=2,k=2"
         i = 2; j = 2; k = 2
         do IR=2,NR-1 ! omit halo elements / nodes
             write(143,'(2f16.4)') ROD_RADIUS_NODES(IR),   tf(i,j,k,IR)        
             write(143,'(2f16.4)') ROD_RADIUS_NODES(IR+1), tf(i,j,k,IR)        
         enddo
         close(143)

      end subroutine write_tf

!         call write_water_temp()
!         open(143, file='water-temp.dat')
!         write(143,*) "# x y z nPhase value"
!         do i=1,NX
!         do j=1,NY
!         do k=1,NZ
!         do iphase=1,NPHASE
!             write(143,'(4i8,f16.4)') i, j, k, iphase, tw(i,j,k,iphase)        
!         enddo
!         enddo
!         enddo
!         enddo
!         close(143)
!

