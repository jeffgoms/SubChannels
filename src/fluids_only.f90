      program fluids_only
      
         implicit none
       INTEGER, PARAMETER :: NX=3,NY=3,NZ=10,NR=6,NPHASE=1!,NPHASE=2
!       INTEGER, PARAMETER :: NX=3,NY=3,NZ=62,NR=32,NPHASE=1
       INTEGER, PARAMETER :: NX_WATER=NX+1,NY_WATER=NY+1,NZ_WATER=NZ
         INTEGER :: TF_MAX_ITS, TW_MAX_ITS, NSTEP_ITS
         REAL :: DX(3),ROD_RADIUS_NODES(NX,NY,NZ,NR+1),ROD_RADIUS,DT
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
!         REAL :: PSI(NX,NY,NZ,NG),PSIOLD(NX,NY,NZ,NG), S_NEUTRONS(NX,NY,NZ,NG)
!         REAL :: SPEED(NG),NEUTRON_DIFF(NX,NY,NZ,NG)
!         REAL :: C(NX,NY,NZ,NDELAY), COLD(NX,NY,NZ,NDELAY)
!         REAL :: BETA_EFF,LAMBDA(NDELAY),BETA(NDELAY)
!         REAL :: SIGMA_MATRIX(NX,NY,NZ,NG,NG),SIGMA_F(NX,NY,NZ,NG), XI_P(NX,NY,NZ,NG),XI_D(NX,NY,NZ,NG),NU_F(NX,NY,NZ,NG)
!         INTEGER :: MATERIAL(NX,NY,NZ)
!         REAL :: INTERP_TEMPERATURE(MAX_NX_TI,NDIM_INTERP )
!         REAL :: NEUTRON_DIFF_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2)
!         REAL ::  SIGMA_A_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4), SIGMA_S_TEMPERATURE(NG,NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
!         REAL :: SIGMA_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
!         REAL :: SIGMA_F_MATRIX_ZERO(NX,NY,NZ,NG,NG), SIGMA_F_MATRIX(NX,NY,NZ,NG,NG)
!         REAL :: XI_P_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4),XI_D_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
!         REAL :: NU_F_TEMPERATURE(NG,NMATERIAL,NX_T1,NX_T2, NX_T3,NX_T4)
         REAL :: VOLFRA_ERROR_TOLER, VEL_ERROR_TOLER, NEU_ERROR_TOLER, NEU_ERROR_TOLER_IN_G
         INTEGER :: VOLFRA_MAX_ITS, VEL_MAX_ITS, NEU_MAX_ITS, NEU_MAX_ITS_IN_G
         INTEGER :: VOLFRA_RELAX, VEL_RELAX, NEU_RELAX
         INTEGER N_RADIAL_ROD_MATERIAL
!         INTEGER, DIMENSION(:,:,:,:), allocatable :: RADIAL_ROD_MATERIAL
!         REAL, DIMENSION(:), allocatable :: MATERIAL_DEN_CP_FUEL, MATERIAL_DIFF_FUEL
!! Local variables WITHIN THIS PROGRAMME...
!         REAL :: EIG, EIG_PREV, EIG_ERROR_TOLER
!         INTEGER :: ITIME, I,J,K,IR, EIG_IT, MAX_EIG_ITS, NTIME
!         REAL :: SIGMA_F_KEEP(NX,NY,NZ,NG), PSI2(NX,NY,NZ,NG)
         REAL, DIMENSION(:,:, :,:), allocatable :: S_FUEL_RODS
         INTEGER ITIME, NTIME, I, j, k, iphase
!        DO ITS=1,NSTEP_ITS  !  ITERATE WITHIN A TIME STEP IF NEEDED. 

! work out what the power source is
        ALLOCATE(S_FUEL_RODS(NX,NY,NZ,NR)); S_FUEL_RODS=0.0


! UNITS IN cm and g
      
        DX(1)=0.02 ! SIZE OF ROD CENTRED CELL IN X DIRECTION
        DX(2)=DX(1)
        DX(3)= 0.05! SIZE OF ROD CENTRED CELL IN Z-DIRECTION
        VOLFRA=1.0;VOLFRAOLD=1.0
!        VEL=100.0; VELOLD=VEL ! Vel of water in m/s
        VEL=1.0; VELOLD=VEL ! Vel of water in m/s
        DEN_WATER=1.0e+3; den_water_old=den_water
        DEN_WATER_BCS=1.0e+3
        DEN_CP_FUEL_RODS=10.0e+3
        DIFF_FUEL_RODS=1.0e+2 ! conductitivity of fuel rods
!        DIFF_FUEL_RODS=1.0 ! conductitivity of fuel rods
!        S_FUEL_RODS=21000000.0 ! SOURCE OF ENERGY W/CM3
        S_FUEL_RODS=21.0 ! SOURCE OF ENERGY W/M3
!        S_FUEL_RODS=0.0 ! SOURCE OF ENERGY W/M3
         
        TF_RELAX=3 ! SWITCH ON TRI-DIAGONAL SOLVE
        TF_MAX_ITS=100
        TF_ERROR_TOLER=1.E-7
        TW_RELAX=1
        TW_MAX_ITS=100
        TW_ERROR_TOLER=1.E-7

         TFOLD = 0.0
         TF = TFOLD
         TWOLD = 0.0
         TW = TWOLD

         p=0.0
         P_BCS=0.0
         TW_BCS=200.0 ! SET INLET TEMP BC
         !ROD_RADIUS_NODES ! NODE POSITIONS IN THE R-DIRECTION
         ROD_RADIUS = 0.007
         ROD_RADIUS_NODES(:,:,:,1:2)=0.0
         DO I=3,NR
           ROD_RADIUS_NODES(:,:,:,I) = REAL(I-2)*ROD_RADIUS/REAL(NR-2)
         ENDDO 
         ROD_RADIUS_NODES(:,:,:,NR+1)=ROD_RADIUS_NODES(:,:,:,NR)

         DT=1.0 ! TIME STEP SIZE SECONDS
         NTIME = 50


         
         print *, 'Running...'

        DO ITIME=1,NTIME
           print *, 'ITIME=', ITIME
! CALCULATE TF, TW (FUEL ROD AND WATER TEMPS) GIVEN OLD TIME STEP VALUES OF THEM, TFOLD, TWOLD

        call TEMP_CALC(TF,TFOLD, TW,TWOLD, VOLFRA,VOLFRAOLD,DEN_WATER,DEN_WATER_OLD,P,VEL,VELOLD, &
                             TW_BCS, DEN_WATER_BCS, VEL_BCS, VOLFRA_BCS, P_BCS, P_BCS_ON, &
                             DX,ROD_RADIUS_NODES,DT,NX,NY,NZ,NX_WATER,NY_WATER,NZ_WATER,NR,NPHASE, DEN_CP_FUEL_RODS, DIFF_FUEL_RODS, S_FUEL_RODS, &
                             TF_RELAX,TF_MAX_ITS,TF_ERROR_TOLER, &
                             TW_RELAX,TW_MAX_ITS,TW_ERROR_TOLER, &
                             VOLFRA_RELAX,VOLFRA_MAX_ITS,VOLFRA_ERROR_TOLER, &
                             VEL_RELAX,VEL_MAX_ITS,VEL_ERROR_TOLER)
         END DO


         print *, "outputting results"



         ! some routines which output cell number and temperatures ... 
         ! writes temp of water in the horizontal plane at the bottom
         call write_tw_btm(NX_WATER, NY_WATER, NZ_WATER, NPHASE, TW) 
         ! writes temp of water in a line in the z-direction
         call write_tw_height(NX_WATER, NY_WATER, NZ_WATER, NPHASE, TW) 
         ! writes temp of water in the horizontal plane at the top
         call write_tw_top(NX_WATER, NY_WATER, NZ, NPHASE, TW) ! should merge this with the previous subroutine ;-)
         ! writes temp of fuel for cell 2,2,2
         call write_tf(NX, NY, NZ, NR, TF, ROD_RADIUS_NODES)



         print *, 'finished'

         end program fluids_only
