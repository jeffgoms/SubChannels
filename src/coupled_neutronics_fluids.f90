       PROGRAM COUPLED_NEUTRONICS_FLUIDS
       IMPLICIT NONE
! Variables to be defined and passed to other 
       INTEGER, PARAMETER :: NX=10,NY=10,NZ=10,NR=5,NPHASE=1, NG=2, NDELAY=6
       INTEGER, PARAMETER :: NX_WATER=NX+1,NY_WATER=NY+1,NZ_WATER=NZ
       INTEGER, PARAMETER :: NDIM_INTERP = 4 ! This is the dimension of the interpolation...
       INTEGER, PARAMETER :: NMATERIAL=1, NX_T1=5, NX_T2=5, NX_T3=1, NX_T4=1, MAX_NX_TI=MAX(NX_T1, NX_T2, NX_T3, NX_T4)
         INTEGER :: TF_MAX_ITS, TW_MAX_ITS, NSTEP_ITS
         REAL :: DX(3),ROD_RADIUS_NODES(NX,NY,NZ,NR+1),DT
         REAL :: TF(NX,NY,NZ,NR),TFOLD(NX,NY,NZ,NR), TW(NX_WATER,NY_WATER,NZ_WATER,NPHASE),TWOLD(NX_WATER,NY_WATER,NZ_WATER,NPHASE)
         REAL :: VOLFRA(NX_WATER,NY_WATER,NZ_WATER,NPHASE),VOLFRAOLD(NX_WATER,NY_WATER,NZ_WATER,NPHASE)
         REAL :: DEN_WATER(NX_WATER,NY_WATER,NZ_WATER,NPHASE),DEN_WATER_OLD(NX_WATER,NY_WATER,NZ_WATER,NPHASE)
         REAL :: P(NX_WATER,NY_WATER,NZ_WATER)
         REAL :: VEL(NX_WATER,NY_WATER,NZ_WATER+1,NPHASE),VELOLD(NX_WATER,NY_WATER,NZ_WATER+1,NPHASE)
         REAL :: TW_BCS(NX_WATER,NY_WATER,NZ_WATER,NPHASE), DEN_WATER_BCS(NX_WATER,NY_WATER,NZ_WATER,NPHASE) 
         REAL :: VEL_BCS(NX_WATER,NY_WATER,NZ_WATER+1,NPHASE), VOLFRA_BCS(NX_WATER,NY_WATER,NZ_WATER,NPHASE), P_BCS(NX_WATER,NY_WATER,NZ_WATER,NPHASE)
         INTEGER :: P_BCS_ON(NX_WATER,NY_WATER,NZ_WATER,NPHASE)
         REAL :: DEN_CP_FUEL_RODS(NX,NY,NZ,NR), DIFF_FUEL_RODS(NX,NY,NZ,NR)
         REAL :: TF_ERROR_TOLER,TW_ERROR_TOLER
         INTEGER :: TF_RELAX,TW_RELAX
         REAL :: PSI(NX,NY,NZ,NG),PSIOLD(NX,NY,NZ,NG), S_NEUTRONS(NX,NY,NZ,NG)
         REAL :: SPEED(NG),NEUTRON_DIFF(NX,NY,NZ,NG)
         REAL :: C(NX,NY,NZ,NDELAY), COLD(NX,NY,NZ,NDELAY)
         REAL :: BETA_EFF,LAMBDA(NDELAY),BETA(NDELAY)
         REAL :: SIGMA_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F(NX,NY,NZ,NG), XI_P(NX,NY,NZ,NG),XI_D(NX,NY,NZ,NG),NU_F(NX,NY,NZ,NG)
         INTEGER :: MATERIAL(NX,NY,NZ)
         REAL :: INTERP_TEMPERATURE(MAX_NX_TI,NDIM_INTERP )
         REAL :: NEUTRON_DIFF_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2)
         REAL ::  SIGMA_A_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4), SIGMA_S_TEMPERATURE(NG,NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL :: SIGMA_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL :: SIGMA_F_MATRIX_ZERO(NX,NY,NZ,NG,NG), SIGMA_F_MATRIX(NX,NY,NZ,NG,NG)
         REAL :: XI_P_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),XI_D_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL :: NU_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL :: VOLFRA_ERROR_TOLER, VEL_ERROR_TOLER, NEU_ERROR_TOLER, NEU_ERROR_TOLER_IN_G
         INTEGER :: VOLFRA_MAX_ITS, VEL_MAX_ITS, NEU_MAX_ITS, NEU_MAX_ITS_IN_G
         INTEGER :: VOLFRA_RELAX, VEL_RELAX, NEU_RELAX
         INTEGER N_RADIAL_ROD_MATERIAL
         INTEGER, DIMENSION(:,:,:,:), allocatable :: RADIAL_ROD_MATERIAL
         REAL, DIMENSION(:), allocatable :: MATERIAL_DEN_CP_FUEL, MATERIAL_DIFF_FUEL
! Local variables WITHIN THIS PROGRAMME...
         REAL :: EIG, EIG_PREV, EIG_ERROR_TOLER
         INTEGER :: ITIME, I,J,K,IR, EIG_IT, MAX_EIG_ITS, NTIME
         REAL :: SIGMA_F_KEEP(NX,NY,NZ,NG), PSI2(NX,NY,NZ,NG)

! ***********************************************************************************************
! *********START DEFINE VARIABLES****************************************************************
! ***********************************************************************************************
         print *, 'coupled neutronics and fluids calculation'
         print *, 'fluids and thermal properties'
! Fluids and thermal properties...
         VEL=300.0; VELOLD=300.0 ! Velocity upwards past the fuel rods - for single phase it is not solved for.
         VOLFRA=1.0;VOLFRAOLD=1.0 ! Water volume fraction - for single phase it is not solved for.
         TF=200.0; TFOLD=200.0; TW=200.0; TWOLD=200.0 ! Fuel and water temp.
         DEN_WATER=1.0; DEN_WATER_OLD=1.0 ! Water density.
!
         DT=1.0 ! Time step size. 
         DX(1)=10.0; DX(2)=10.0 ! 10 cm between fuel rods
         DX(3)=10.0 ! Distance between CV spacing in the vertical.
         ROD_RADIUS_NODES(:,:,:,1)=0.0  ! The control rod is 3cm in radius. 
         ROD_RADIUS_NODES(:,:,:,2)=0.0
         ROD_RADIUS_NODES(:,:,:,3)=1.0
         ROD_RADIUS_NODES(:,:,:,4)=2.0
         ROD_RADIUS_NODES(:,:,:,NR)=3.0
         ROD_RADIUS_NODES(:,:,:,NR+1)=3.0

         P=0.0 ! Bulk pressure - for single phase it is not solved for.
! IMAT = RADIAL_ROD_MATERIAL(I,J,K,IR) - contains the materials for the fuel integer
         !RADIAL_ROD_MATERIAL=1 ! JUST ONE MATERIAL
         N_RADIAL_ROD_MATERIAL=1 ! JUST ONE MATERIAL


         ALLOCATE(RADIAL_ROD_MATERIAL(NX,NY,NZ,NR),MATERIAL_DEN_CP_FUEL(N_RADIAL_ROD_MATERIAL),MATERIAL_DIFF_FUEL(N_RADIAL_ROD_MATERIAL))
         RADIAL_ROD_MATERIAL=1 ! Assume rod material 1 is fissile -needed for distributing the source of heat in fuel
         MATERIAL_DEN_CP_FUEL=1.0
         MATERIAL_DIFF_FUEL=1.0
         DO IR=1,NR
         DO K=1,NZ
         DO J=1,NY
         DO I=1,NX
            DEN_CP_FUEL_RODS(I,J,K,IR) = MATERIAL_DEN_CP_FUEL( RADIAL_ROD_MATERIAL(I,J,K,IR) ) ! density*cp of fuel
            DIFF_FUEL_RODS(I,J,K,IR) = MATERIAL_DIFF_FUEL( RADIAL_ROD_MATERIAL(I,J,K,IR) )! conductivity of fuel
         END DO
         END DO
         END DO
         END DO

         print *, 'neutronics properties'
! Neutronics properties...
! Assume a cell with NEUTRON_DIFF(I,J,K,G) <= 0.0 is outside the calculations domain
! =0 assume a corresponding vacume bc <0 assume a symmetry b.c.
         S_NEUTRONS=0.0 ! Source of neutrons

! fuel rod properties...       
         DEN_CP_FUEL_RODS = 0.0 ! density * C_p for fuel/control rods
         DIFF_FUEL_RODS = 0.0 ! temp diffusion coeff for fuel/control rods

! define bcs: 
! TW_BCS, DEN_WATER_BCS, VOLFRA_BCS, P_BCS contains the bcs of water temp, water density, volume fraction and pressure.
! P_BCS_ON =1 switches on the pressure bc for a particular phase, P_BCS_ON =0 assume that we have a vel b.c.
         TW_BCS = 200.0; DEN_WATER_BCS = 1.0; VOLFRA_BCS = 1.0; P_BCS = 0.0
         P_BCS_ON = 0
         K=1 ! Bottom
         TW_BCS(:,:,K,:) = 200.0; DEN_WATER_BCS(:,:,K,:) = 1.0; VOLFRA_BCS(:,:,K,1) = 1.0; VOLFRA_BCS(:,:,K,2:NPHASE) = 0.0 
         VEL_BCS(:,:,K,:) = 0.0
         P_BCS(:,:,K,:) = 1.E+6
         P_BCS_ON(:,:,K,:)=0
         K=NZ_WATER ! Top
         TW_BCS(:,:,K,:) = 200.0; DEN_WATER_BCS(:,:,K,:) = 1.0; VOLFRA_BCS(:,:,K,1) = 1.0; VOLFRA_BCS(:,:,K,2:NPHASE) = 0.0 
         VEL_BCS(:,:,K+1,:) = 0.0
         P_BCS(:,:,K,:) = 1.E+6
         P_BCS_ON(:,:,K,:)=1

! SOLVER OPTIONS: _WATER
! Rod temp solver:
         TF_ERROR_TOLER = 1.e-7
         TF_MAX_ITS = 10
         TF_RELAX = 0 ! Use tridiagonal solver
! Water temp solver...
         TW_ERROR_TOLER = 1.e-7
         TW_MAX_ITS = 10
         TW_RELAX = 1; IF(NPHASE.GT.1) TW_RELAX = 0 ! Forward Backward Block Gauss-Siedel solver

! VOLFRA soln (volume fraction): 
         VOLFRA_ERROR_TOLER = 1.e-7
         VOLFRA_MAX_ITS = 20
         VOLFRA_RELAX= 1; IF(NPHASE.GT.1) VOLFRA_RELAX = 3 ! Forward Backward Block Gauss-Siedel solver
! VEL solver:
         VEL_ERROR_TOLER = 1.e-7
         VEL_MAX_ITS = 20
         VEL_RELAX= 1; IF(NPHASE.GT.1) VEL_RELAX = 3 ! Forward Backward Block Gauss-Siedel solver

! Neutronics solution:
         NEU_ERROR_TOLER = 1.e-7 ! Outer iteration between groups
         NEU_MAX_ITS = 10 ! Outer iteration between groups
         NEU_ERROR_TOLER_IN_G = 1.e-7 ! Inner iteration within a group
         NEU_MAX_ITS_IN_G = 10 ! Inner iteration within a group
         NEU_RELAX = 0 ! tri-diagonal solver - block diagonal solver
         EIG_ERROR_TOLER = 1.E-7

! ITERATE WITHIN A TIME STEP IF NEEDED...
         NSTEP_ITS=1 
         MAX_EIG_ITS=100
         NTIME=0 ! =0 FOR eigenvalue problem.
! ***********************************************************************************************
! *********END DEFINE VARIABLES****************************************************************
! ***********************************************************************************************



       print *, 'running...'
       IF(NTIME.GT.0) THEN
          DO ITIME=1,NTIME
             print *, 'time stepping...',itime,' of ',ntime
             CALL NEUTRONICS_STEP(TF,TFOLD, TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P,VEL,VELOLD, &
                             TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS, P_BCS_ON, &
                             DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NX_WATER,NY_WATER,NZ_WATER,NG,NR,NPHASE, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, RADIAL_ROD_MATERIAL,&
                             TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, &
                             TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, &
                             PSI,PSIOLD,C,COLD,S_NEUTRONS,SPEED,NEUTRON_DIFF,SIGMA_MATRIX, NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER,NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G, &
                             NDELAY,BETA_EFF,BETA,LAMBDA,XI_D,XI_P,SIGMA_F_MATRIX,SIGMA_F,NU_F, &
                             NSTEP_ITS, & 
                             INTERP_TEMPERATURE, NEUTRON_DIFF_TEMPERATURE,SIGMA_A_TEMPERATURE, SIGMA_S_TEMPERATURE, SIGMA_F_TEMPERATURE, &
                             XI_P_TEMPERATURE,XI_D_TEMPERATURE,NU_F_TEMPERATURE, MATERIAL,NMATERIAL, NX_T1,NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP, &
                             VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                             VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER, NTIME)
          END DO
       ELSE ! Eigen problem
! INTERPOLATE X-SECTIONS IN TEMPERATURE...
           print *, 'eigen problem'           
           CALL INTERPOLATE_XSECTIONS_4D(TF,TW,VOLFRA,DEN_WATER,NEUTRON_DIFF,SIGMA_MATRIX,SIGMA_F_MATRIX,SIGMA_F, XI_P,XI_D,NU_F, ROD_RADIUS_NODES, &
                            NX,NY,NZ,NX_WATER,NY_WATER,NZ_WATER,NG,NR,NDELAY,NPHASE, &
                            INTERP_TEMPERATURE, NEUTRON_DIFF_TEMPERATURE,SIGMA_A_TEMPERATURE, SIGMA_S_TEMPERATURE, SIGMA_F_TEMPERATURE, &
                            XI_P_TEMPERATURE,XI_D_TEMPERATURE,NU_F_TEMPERATURE, MATERIAL,NMATERIAL, NX_T1,NX_T2, NX_T3,NX_T4, MAX_NX_TI, NDIM_INTERP, &
                            BETA_EFF, NTIME )

          EIG_PREV=0.0
          PSI=1.0
          DO EIG_IT=1,MAX_EIG_ITS

             PSI2=PSI
             S_NEUTRONS=0.0
! MATRIX VECTOR FOR S_NEUTRONS=FISSION_MATRIX*PSI
             CALL FISSION_MATRIX_VEC(PSI,SIGMA_F_KEEP,S_NEUTRONS, NX,NY,NZ,NG)

             SIGMA_F_MATRIX_ZERO=0.0  
! This sub calculates the scalar flux PSI for each energy group. 
             CALL NEUTRONICS_DIFFUSION(PSI,PSIOLD,S_NEUTRONS,SPEED,NEUTRON_DIFF,SIGMA_MATRIX,DT,DX, NX,NY,NZ,NG, &
                                     NEU_RELAX,NEU_MAX_ITS,NEU_ERROR_TOLER,NEU_MAX_ITS_IN_G,NEU_ERROR_TOLER_IN_G, &
                                     NDELAY,BETA_EFF,BETA,LAMBDA,XI_D,XI_P,SIGMA_F_MATRIX_ZERO,SIGMA_F,C) 

             EIG = SUM(PSI*PSI2)/SUM(PSI2**2)
             IF(ABS(EIG-EIG_PREV).LT.EIG_ERROR_TOLER) CYCLE
             EIG_PREV=EIG

          END DO ! DO EIG_IT=1,NEIG_ITS

          PSI = PSI/SUM(PSI2**2) ! Normalize

       ENDIF
        
        print *, 'finishing coupled neutronics fluids calculation' 
        RETURN
        END PROGRAM COUPLED_NEUTRONICS_FLUIDS ! end of main program for coupled neutronics-fluids. 

