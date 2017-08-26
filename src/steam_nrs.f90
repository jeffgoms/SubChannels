!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details. 
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
!#include "fdebug.h"

module steam_nrs_module

!  use FLDebug
!  use AllSorts

  contains

SUBROUTINE STEAM_NRS(P,T,DENSTF)
 
  implicit none
  real t,p,denstf

  call outvolume(t,p,denstf)
  
return
end subroutine

subroutine outvolume(t,p,rho)
  !
  !*******************************************************************************
  !
  !! outvolume tests VOLUME.
  !
  !
  real dpdr
  real dvdr
  real dvdt
  integer i
  integer j
  real pmpa
  real rho
  real rhostart
  real tk
  real v

  real p, t
  real pp, tt, tsatur
  real rhol, rhov, xv

!!! Units conversion:
!!!   1. Pressure:
!!!       Input: g/(cm.s^2) --> Output: MPa
!!!   2. Temperature:
!!!       Input: degree C   --> Output: K

  tt = t + 273.15
  pp = p * 1.e-7
  rhostart = 1.9 ; rhol = rhostart ; rhov = rhostart/1000. ! Initialisation of density

  call dense(pp,tt,rhostart, rho, dpdr )
  call calc_volume(tt,rho,v,dvdt,dvdr)
  call tsat( pp, tsatur, rhol, rhov )
  print*, 'Tsat:', tsatur
  !print*, 'rhol, rhov:', 1./rhol, 1./rhov

!!! Calculating the steam quality
   xv = min( 1., max( 0., (rho - rhol)/(rhov-rhol) ) )
   print*, xv, rho, rhol, rhov
   print*, xv*rhov, (1.-xv)*rhol


  return
end subroutine outvolume



subroutine base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )
!
!*******************************************************************************
!
!! BASE calculates quantities associated with the base Helmholtz function.
!
!
!  Discussion:
!
!    The equation for the base Helmholtz function AB(T,RHO) is:
!
!      AB(T,RHO) = R * T * (
!        - ln ( 1 - y ) 
!        - ( beta - 1 ) / ( 1 - y ) 
!        + ( alpha + beta + 1 ) / ( 2 * ( 1 - y )**2 )
!        + 4 * y * ( ( Bbar / b ) - gamma ) 
!        - 0.5 * ( alpha - beta + 3 ) 
!        + ln ( RHO * R * T / P0 ) )
!                                                      (Equation 2)
!   where 
!
!     y = b * rho / 4, 
!     alpha = 11,
!     beta = 133/3, 
!     gamma = 7/2, 
!     P0 = 0.101325 MegaPascals = 1 atm
!
!   and
!
!     b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!     Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!   where 
!
!     T0=647.073 K and the coefficients b(j) and B(j) are
!    
!     j    b(j)                         B(j)
!    --    -----------                  ----------
!     0    0.7478629                    1.1278334
!     1   -0.3540782                   -0.5944001
!     2    0                           -5.010996
!     3    0.007159876                  0
!     4    0                            0.63684256
!     5   -0.003528426                  0
!
!  For the derived quantities, the following relations are used:
!
!    Pressure:                  PB      = RHO**2 * dAB/dRHO
!    Density derivative:        DPDRB   = 2*PB/RHO + RHO**2 * d2AB/dRHO2
!    Temperature derivative:    DPDTB   = RHO**2 * d2AB/(dRHO dT)
!    Specific entropy:          SB      = ( UB - AB ) / T
!    Specific internal energy:  UB      = AB + T * SB
!    Specific enthalpy:         HB      = UB + PB / RHO
!    Specific heat capacity
!      at constant volume:      CVB     = - T * d2AB/dT2
!    Specific Gibbs function:   GB      = AB + PB / RHO
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the density, in G/CM3.
!
!    Output, real AB, the base value of the Helmholtz function,
!    in KJ/KG.
!
!    Output, real CVB, the base value of the isochoric (constant 
!    volume) heat capacity, in KJ/(KG degrees Kelvin).
!
!    Output, real DPDRB, the base value of the partial
!    derivative dP(T,RHO)/dRHO, with T held fixed, in (MegaPascals CM3)/G.
!
!    Output, real DPDTB, the base value of the partial
!    derivative dP(T,RHO)/dT, with RHO held fixed, in MegaPascals/degrees Kelvin.
!
!    Output, real GB, the base value of the Gibbs free energy,
!    in KJ/KG.
!
!    Output, real HB, the base value of enthalpy, in KJ/KG.
!
!    Output, real PB, the base pressure, in MegaPascals.
!
!    Output, real SB, the base value of entropy, 
!    in KJ/(KG degrees Kelvin).
!
!    Output, real UB, the base value of internal energy, 
!    in KJ/KG.
!
  real, parameter :: ALPHA = 11.0D+00
  real, parameter :: BETA = 44.333333333333D+00
  real, parameter :: GAMMA = 3.5D+00
  real, parameter :: PZERO = 0.101325D+00
!
  real ab
  real b1
  real b1t
  real b1tt
  real b2
  real b2t
  real b2tt
  real cvb
  real dpdrb
  real dpdtb
  real dz
  real dz0
  real gb
  real hb
  real pb
  real rho
  real sb
  real t
  real ub
  real x
  real y
  real z
  real z0
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0 ) then
    write(*, * ) ' '
    write(*, * ) 'BASE - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'BASE - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute auxilliary quantities for Equation 2.
!
  call bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )

  y = 0.25D+00 * b1 * rho

  x = 1.0D+00 - y
!
!  Evaluate Equation 2.
!
  ab =   - log ( 1.0D+00 - y ) &
         - ( BETA - 1.0D+00 ) / ( 1.0D+00 - y ) &
         + ( ALPHA + BETA + 1.0D+00 ) / ( 2.0D+00 * ( 1.0D+00 - y )**2 ) &
         + 4.0D+00 * y * ( ( b2 / b1 ) - GAMMA ) &
         - 0.5D+00 * ( ALPHA - BETA + 3.0D+00 ) &
         + log ( rho * calc_gascon() * t / PZERO )
!
!  Determine quantities defined in terms of AB.
!
  pb = ( 1.0D+00 + ALPHA * y + BETA * y**2 ) / ( 1.0D+00 - y )**3 &
    + 4.0D+00 * y * ( b2 / b1 - GAMMA )

  z0 = ( 1.0D+00 + ALPHA * y + BETA * y * y ) / x**3

  z = z0 + 4.0D+00 * y * ( b2 / b1 - GAMMA )

  dz0 = ( ALPHA + 2.0D+00 * BETA * y ) / ( 1.0D+00 - y )**3 &
    + 3.0D+00 * ( 1.0D+00 + ALPHA * y + BETA * y**2 ) / ( 1.0D+00 - y )**4

  dz = dz0 + 4.0D+00 * ( b2 / b1 - GAMMA )

  gb = ab + pb

  ub = - t * b1t * ( pb - 1.0D+00 - rho * b2 ) / b1 - rho * t * b2t

  hb = pb + ub

  cvb = 2.0D+00 * ub + ( z0 - 1.0D+00 ) &
    * ( ( t * b1t / b1 )**2 - t**2 * b1tt / b1 ) &
    - rho * t**2 * ( b2tt - GAMMA * b1tt ) - ( t * b1t / b1 )**2 * y * dz0

  dpdtb = pb / t + rho * ( 0.25D+00 * ( dz0 + 4.0D+00 * ( b2 / b1 - GAMMA ) ) &
    * b1t + b2t - b2 / b1 * b1t )

  sb = ub - ab

  dpdrb = pb + y * ( dz0 + 4.0D+00 * ( b2 / b1 - GAMMA ) )
!
!  Assign dimensions.
!
  ab =    calc_gascon() * t       * ab
  cvb =   calc_gascon()           * cvb
  dpdrb = calc_gascon() * t       * dpdrb
  dpdtb = calc_gascon() * t * rho * dpdtb
  gb =    calc_gascon() * t       * gb
  hb =    calc_gascon() * t       * hb
  pb =    calc_gascon() * t * rho * pb
  sb =    calc_gascon()           * sb
  ub =    calc_gascon() * t       * ub

  return
end subroutine base

subroutine bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )
!
!*******************************************************************************
!
!! BB calculates the B's of equations 3 and 4.  
!
!
!  Discussion:
!
!    Here
!
!      b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!      Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!    where 
!
!      T0=647.073 K and the coefficients b(j) and B(j) are
!    
!      j    b(j)                         B(j)
!     --    -----------                  ----------
!      0    0.7478629                    1.1278334
!      1   -0.3540782                   -0.5944001
!      2    0                           -5.010996
!      3    0.007159876                  0
!      4    0                            0.63684256
!      5   -0.003528426                  0
!
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Output, real B1, the coefficient b from equation 3, 
!    in CM3/G.
!
!    Output, real B2, the coefficient Bbar from equation 4, 
!    in CM3/G.
!
!    Output, real B1T, the derivative dB1/dT, 
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, real B2T, the derivative dB2/dT, 
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, real B1TT, the second derivative of B1 with 
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
!    Output, real B2TT, the second derivative of B2 with 
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
  real, parameter :: TREF = 647.073D+00
!
  real b1
  real b1t
  real b1tt
  real b2
  real b2t
  real b2tt
  real, parameter, dimension ( 10 ) :: bp = (/ &
    0.7478629D+00,   -0.3540782D+00,    0.0D+00,           0.0D+00, &
    0.007159876D+00,  0.0D+00,         -0.003528426D+00,   0.0D+00, &
    0.0D+00,          0.0D+00 /)
  real, parameter, dimension ( 10 ) :: bq = (/ &
    1.1278334D+00,    0.0D+00,         -0.5944001D+00,   -5.010996D+00, &
    0.0D+00,          0.63684256D+00,   0.0D+00,          0.0D+00, &
    0.0D+00,          0.0D+00 /)
  integer i
  real t
  real v(10)
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'BB - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Set V(I) = ( TREF / T )**(I-1).
!
  v(1) = 1.0D+00
  do i = 2, 10
    v(i) = v(i-1) * TREF / t
  end do
!
!  Set B1, B1T, B1TT.
!
  b1 = bp(1) + bp(2) * log ( 1.0D+00 / v(2) )
  b1t = bp(2) * v(2) / TREF
  b1tt = 0.0D+00
  do i = 3, 10
    b1 = b1 + bp(i) * v(i-1)
    b1t = b1t - dble ( i - 2 ) * bp(i) * v(i-1) / t
    b1tt = b1tt + bp(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b1tt = b1tt -  ( b1t / t )
!
!  Set B2, B2T, B2TT.
!
  b2 = bq(1)
  b2t = 0.0D+00
  b2tt = 0.0D+00
  do i = 3, 10
    b2 = b2 + bq(i) * v(i-1)
    b2t = b2t - dble ( i - 2 ) * bq(i) * v(i-1) / t
    b2tt = b2tt + bq(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b2tt = b2tt - ( b2t / t )

  return
end subroutine bb

subroutine corr ( t, p, rhol, rhov, delg )
!
!*******************************************************************************
!
!! CORR evaluates an adjustment to the Gibbs function.
! 
!
!  Discussion:
!
!    CORR is given T and P at or near the vapor pressure and evaluates
!    the corresponding liquid and vapor densities, and the residual
!    function DELG = (GL-GV)/(R*T) where GL and GV are the Gibbs functions
!    for the liquid and vapor phases, respectively.
!
!    These quantities are used to calculate a correction to the vapor 
!    pressure or the vapor temperature.  
!
!    The states corresponding to the coexisting phases of liquid
!    and vapor for the temperature range from the triple point
!    to within 0.5 C of the critical point 0.01 <= t <= tk-0.5 C
!    have been determined in exact accord with the Gibbs condition
!    of phase equilibrium: DELG = G(g)-G(l) = 0, P, t constant,
!    where G(g) and G(l) are the values of the Gibbs function
!    for saturated gas and liquid respectively.
!   
!    For the region (tk-t)<=0.5 C, an exact solution for the
!    Helmholtz function yields values of density for the saturated
!    liquid that are shifted to lower values.  Also, the isotherms
!    in the pressure-density plane and the Gibbs function-density
!    plane are nearly flat, so that it is difficult to obtain
!    solutions.  As an alternative to exact solution, the power
!    law equation is used to define states:
!  
!      rho(gas) = 0.322 - 0.657 * (1 - T/647.126)**0.325   (g/cm3).
!      rho(liq) = 0.322 + 0.657 * (1 - T/647.126)**0.325   (g/cm3).
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the vapor temperature, in degrees Kelvin.
!
!    Input/output, real P, the vapor pressure, in MegaPascals.
!    On output, if 646.3< = T <= TCRIT, P is replaced
!    by a value consistent with the computed RHOL and RHOV.
!    For TCRIT<T, P is set to 22.055.
!
!    Input/output, real RHOL, the liquid density, in G/CM3.
!    On input, if RHOL is greater than 0, it is used as an initial
!    estimate for the iteration.
!
!    Input/output, real RHOV, the vapor density, in G/CM3.
!    On input, if RHOV is greater than 0, it is used as an initial
!    estimate for the iteration.
!
!    Output, real DELG, the residual function (GL-GV)/(R*T),
!    where GL is the liquid Gibbs function, GV the vapor Gibbs function, 
!    dimensionless.  If T > 646.3, DELG is 0.
!
  real, parameter :: PCRIT = 22.055D+00
  real, parameter :: TCRIT = 647.1260000001D+00
!
  real a
  real ab
  real ar
  real cd
  real cjth
  real cjtt
  real cp
  real cv
  real cvb
  real cvr
  logical, parameter :: LOCAL_DEBUG = .false.
  real delg
  real dpdr
  real dpdrb
  real dpdrr
  real dpdt
  real dpdtb
  real dpdtr
  real g
  real gb
  real gl
  real gr
  real gv
  real h
  real hb
  real hr
  real p
  real pb
  real pr
  real psave
  real rho
  real rhol
  real, parameter :: RHO_MIN = 1.0D-08
  real rhov
  real rhostart
  real s
  real sb
  real sr
  real t
  real tau
  real u
  real ub
  real ur
!
  psave = p
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'CORR - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'CORR - Fatal error!'
    write(*, * ) '  Input pressure P must be positive.'
    write(*, * ) '  Input value was P = ', p
    stop
  end if

  if ( t <= 646.3D+00 ) then

    if ( rhol <= 0.0D+00 ) then
      rhostart = 1.11D+00 - 0.0004D+00 * t
    else
      rhostart = rhol
    end if

    call dense ( p, t, rhostart, rho, dpdr )

    call therm ( t, rho, a, cjth, cjtt, cd, cv, dpdr, dpdt, g, h, p, s, u )

    p = psave
    rhol = rho
    gl = g

    if ( rhov <= 0.0D+00 ) then
      rhostart = p / ( calc_gascon() * t )
    else
      rhostart = rhov
    end if

    call dense ( p, t, rhostart, rho, dpdr )

    rho = max ( rho, RHO_MIN )

    call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

    rhov = rho
    gv = g
    delg = ( gl - gv ) / ( calc_gascon() * t )

    p = psave

    if ( LOCAL_DEBUG ) then
      write(*, * ) 'CORR - RHOL = ', rhol, ' RHOV = ', rhov
    end if

  else if ( t <= TCRIT ) then

    if ( LOCAL_DEBUG ) then
      write(*, * ) 'CORR - Twilight zone'
    end if

    delg = 0.0D+00
    tau = 0.657128D+00 * ( 1.0D+00 - t / TCRIT )**0.325D+00
    rhol = 0.322D+00 + tau
    rhov = 0.322D+00 - tau
    rho = rhov
    call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )
    call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )
    p = pb + pr

  else

    if ( LOCAL_DEBUG ) then
      write(*, * ) ' '
      write(*, * ) 'CORR - Weirdo zone'
    end if

    rhol = 0.322D+00
    rhov = 0.322D+00
    p = PCRIT
    delg = 0.0D+00

  end if
 
  return
end subroutine corr

subroutine dense ( p, t, rhostart, rho, dpdr )
!
!*******************************************************************************
!
!! DENSE computes the density for a given pressure and temperature.
!
!
!  Modified:
!
!    19 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real P, the pressure, in MegaPascals.
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHOSTART, an initial guess for the density, 
!    in G/CM3, which also signals whether a vapor or liquid
!    calculation is to be done.  If DPDR is computed negative, then for 
!    0.2967 <= RHOSTART, liquid is assumed, otherwise gas.
!
!    Output, real RHO, the density for the given
!    pressure and temperature, in G/CM3.
!
!    Output, real DPDR, the partial derivative
!    dP(T,RHO)/dRHO, with T held fixed, in (MegaPascals CM3)/G.
!
  integer, parameter :: ITMAX = 50
  real, parameter :: RHO_MAX = 1.9D+00
  real, parameter :: RHO_MIN = 1.0D-08
!
  real ab
  real ar
  real cvb
  real cvr
  real dp
  real dpdr
  real dpdrb
  real dpdrr
  real dpdtb
  real dpdtr
  real dpdx
  real ERRTOL
  real gb
  real gr
  real hb
  real hr
  integer it
  real p
  real pb
  real pp
  real pr
  real rho
  real rhostart
  real sb
  real sr
  real t
  real ub
  real ur
  real x
!
  ERRTOL = sqrt ( epsilon ( ERRTOL ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'DENSE - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'DENSE - Fatal error!'
    write(*, * ) '  Input pressure P must be positive.'
    write(*, * ) '  Input value was P = ', p
    stop
  end if

  rho = rhostart
  rho = max ( rho, RHO_MIN )
  rho = min ( rho, RHO_MAX )

  do it = 1, ITMAX

    call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

    call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

    pp = pb + pr
    dpdr = dpdrb + dpdrr
!
!  Check for negative DP/DRho, which characterizes the two-phase region.
!
    if ( dpdr <= 0.0D+00 ) then

      if ( rhostart >= 0.2967D+00 ) then
        rho = rho * 1.02D+00
      else
        rho = rho * 0.98D+00
      end if

      if ( it <= 10 ) then
        cycle
      end if

    end if

    dpdx = 1.1D+00 * dpdr
    dpdx = max ( dpdx, 0.01 )

    dp = abs ( 1.0D+00 - pp / p )

    if ( dp < ERRTOL .or. &
       ( rho > 0.3D+00 .and. dp < ERRTOL ) .or. &
       ( rho > 0.7D+00 .and. dp < 10.0D+00 * ERRTOL ) ) then
      return
    end if

    x = ( p - pp ) / dpdx
    if ( abs ( x ) > 0.1D+00 ) then
      x = x * 0.1D+00 / abs ( x )
    end if

    rho = rho + x

    rho = max ( rho, RHO_MIN )
    rho = min ( rho, RHO_MAX )

  end do

  write(*, * ) ' '
  write(*, * ) 'DENSE - Warning!'
  write(*, * ) '  The iteration did not converge.'
  write(*, * ) '  Number of iterations was ', ITMAX
  write(*, * ) '  Last iterate was ', rho

  return
end subroutine dense

function calc_gascon ( )
!
!*******************************************************************************
!
!! GASCON returns the value of the specific gas constant.
!
!
!  Note:
!
!    The specific gas constant R is related to the universal gas
!    constant R-bar = 8.31441 J/(mol degrees Kelvin) by the molar mass 
!    M = 18.0152 g/mol:
!
!      R = R-bar / M.
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Output, real GASCON, the value of the specific gas 
!    constant, in J/(g degrees Kelvin).
!
  real calc_gascon
!
  calc_gascon = 0.461522D+00
 
  return
end function calc_gascon

subroutine ideal ( t, ai, cpi, cvi, gi, hi, si, ui )
!
!*******************************************************************************
!
!! IDEAL computes ideal gas thermodynamic properties of water.
!
!
!  Discussion:
!
!    Values for thermodynamic properties of water in the ideal
!    gas state were reported by Woolley.  The formula for the ideal gas
!    term of the Helmholtz function approximates a term by term summation of 
!    contributions from each of the rotation and vibration states.  
!    The formula is
!   
!    A(ideal)(T) = -R*T * (1 + ((C(1)/Tr)+C(2))*ln(Tr)
!      + sum(i=3 to 18) C(i)*Tr**(i-6)                         (Equation 6)
!   
!    where Tr=T/100 K.  The C(i) are tabulated coefficients.  Equation
!    6 can be used for temperatures below 3000 K, and is accurate to
!    within the tolerance of the gas constant for 50<=T<=2000 K.
!     
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Output, real AI, the Helmholtz function, in KJ/KG.
!
!    Output, real CPI, the heat capacity at constant pressure,
!    in KJ/(KG degrees Kelvin).
!
!    Output, real CVI, the heat capacity at constant volume,
!    in KJ/(KG degrees Kelvin).
!
!    Output, real GI, the Gibbs free energy, in KJ/KG.
!
!    Output, real HI, the enthalpy, in KJ/KG.
!
!    Output, real SI, the entropy, in KJ/(KG degrees Kelvin).
!
!    Output, real UI, the internal energy, in KJ/KG.
!
  real ai
  real, parameter, dimension ( 18 ) :: c = (/ &
   19.730271018D+00,      20.9662681977D+00,    -0.483429455355D+00, &
    6.05743189245D+00,    22.56023885D+00,      -9.87532442D+00, &
   -4.3135538513D+00,     0.458155781D+00,      -0.047754901883D+00, &
    0.0041238460633D+00, -0.00027929052852D+00,  0.14481695261D-04, &
   -0.56473658748D-06,     0.16200446D-07,      -0.3303822796D-09, &
    0.451916067368D-11,  -0.370734122708D-13,   0.137546068238D-15 /)
  real cpi
  real cvi
  real gi
  real hi
  integer i
  real si
  real t
  real temp
  real tt
  real ui
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'IDEAL - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if

  tt = t / 100.0D+00
  gi = - ( c(1) / tt + c(2) ) * log ( tt )
  hi = c(2) + c(1) * ( 1.0D+00 - log ( tt ) ) / tt
  cpi = c(2) - c(1) / tt

  do i = 3, 18
    temp = c(i) * tt**(i-6)
    gi = gi - temp
    hi = hi + dble ( i - 6 ) * temp
    cpi = cpi + dble ( ( i - 6 ) * ( i - 5 ) ) * temp
  end do

  ai = gi - 1.0D+00
  ui = hi - 1.0D+00
  cvi = cpi - 1.0D+00
  si = hi - gi
!
!  Assign dimensions.
!
  ai =  calc_gascon() * t * ai
  cpi = calc_gascon()     * cpi
  cvi = calc_gascon()     * cvi
  gi =  calc_gascon() * t * gi
  hi =  calc_gascon() * t * hi
  si =  calc_gascon()     * si
  ui =  calc_gascon() * t * ui

  return
end subroutine ideal

subroutine psat_est ( t, p )
!
!*******************************************************************************
!
!! PSAT_EST makes a rough estimate of the vapor pressure.
!
!
!  Discussion:
!
!    The calculation agrees with tabulated data to within
!    0.02% for temperature to within a degree or so of the critical
!    temperature.  The approximate vapor pressure can be refined
!    by imposing the condition that the Gibbs functions of the vapor
!    and liquid phases be equal.
!
!  Modified:
!
!    21 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Output, real P, the vapor pressure, in MegaPascals.
!
  real, parameter :: TREF = 647.25D+00
!
  real, parameter, dimension ( 8 ) :: a = (/ &
    -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,  -148.39348D+00, 48.631602D+00 /)
  real b
  integer i
  real p
  real q
  real t
  real v
  real w
  real z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'PSAT_EST - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if

  if ( t <= 314.0D+00 ) then

    p = 0.1D+00 * exp ( 6.3573118D+00 - 8858.843D+00 / t &
      + 607.56335D+00 * t**( -0.6D+00 ) )

  else

    v = t / TREF
    w = abs ( 1.0D+00 - v )
    b = 0.0D+00
    do i = 1, 8
      z = i
      b = b + a(i) * w**( ( z + 1.0D+00 ) / 2.0D+00 )
    end do

    q = b / v
    p = 22.093D+00 * exp ( q )

  end if

  return
end subroutine psat_est

subroutine resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )
!
!*******************************************************************************
!
!! RESID calculates residual contributions to thermodynamic quantities.
!
!
!  Discussion:
!
!    The residual function consists of 40 terms.  The first 36 are
!    used in a global least squares fit to experimental data.
!    Three terms were added that contribute only in the immediate
!    neighborhood of the critical point (tk-5)<=t<=(tk+5) C
!    0.20<=p<=0.44 g/cm3, and a single term was added for the
!    region of high pressure and low temperature: t<75 C, P>300 MPa.
!    Except in these limited regions, the residual function is
!    given by the first 36 terms.  The equation is
!   
!      A(residual)(rho,T)=
!        sum(i=1 to 36) (g(i)/k(i)) * (T0/T)**(l(i)) (1-exp(-rho))**(k(i))
!        + sum(i=37 to 40) g(i)*delta(i)**(k(i))
!        * exp(-alpha(i)*delta(i)**(k(i)) - beta(i)*tau(i)**2)
!                                                     (Equation 5)
!   
!    where the g(i) are coefficients determined by fits to data,
!    the delta(i) are reduced densities (delta(i)=((rho-rho(i))/rho(i))
!    and tau(i) are reduced temperatures (tau(i)=((T-tau(i))/tau(i))
!    where rho(i) and tau(i) are specified densities and temperatures.
!    The k(i) and l(i) are specified integers.
!   
!  Modified:
!
!    22 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the density, in G/CM3.
!
!    Output, real AR, the residual contribution to the 
!    Helmholtz function, in KJ/KG.
!
!    Output, real CVR, the residual contribution to the 
!    isochoric (constant volume) heat capacity, in KJ/(KG degrees Kelvin).
!
!    Output, real DPDRR, the residual contribution to 
!    the partial derivative dP(T,RHO)/dRHO, with T held fixed, in 
!    (MegaPascals CM3)/G.
!
!    Output, real DPDTR, the residual contribution to 
!    the partial derivative dP(T,RHO)/dT, with RHO held fixed, 
!    in MegaPascals/degrees Kelvin.
!
!    Output, real GR, the residual contribution to the Gibbs 
!    function, in KJ/KG.
!
!    Output, real HR, the residual contribution to the 
!    enthalpy, in KJ/KG.
!
!    Output, real PR, the residual contribution to the pressure, 
!    in MegaPascals.
!
!    Output, real SR, the residual contribution to the entropy, 
!    in KJ/(KG degrees Kelvin).
!
!    Output, real UR, the residual contribution to the 
!    internal energy, in KJ/KG.
!
  real, parameter :: aa = 1.0D+00
  real, parameter :: SREF = 7.6180720166752D+00
  real, parameter :: TZERO = 647.073D+00
  real, parameter :: UREF = - 4328.4549774261D+00
!
  real, parameter, dimension ( 4 ) :: aad = (/ &
    34.0D+00, 40.0D+00, 30.0D+00, 1050.0D+00 /)
  real, parameter, dimension ( 4 ) :: aat = (/ &
    20000.0D+00, 20000.0D+00, 40000.0D+00, 25.0D+00 /)
  real, parameter, dimension ( 4 ) :: adz = (/ &
    0.319D+00, 0.319D+00, 0.319D+00, 1.55D+00 /)
  real ar
  real att
  real, parameter, dimension ( 4 ) ::  atz = (/ &
    640.0D+00, 640.0D+00, 641.6D+00, 270.0D+00 /)
  real cvr
  real dadt
  real ddz
  real del
  real dex
  real dfdt
  real dpdrr
  real dpdtr
  real e
  real ERRTOL
  real ex0
  real ex1
  real ex2
  real fct
  real, parameter, dimension ( 40 ) :: g = (/ &
    -530.62968529023D+00,  0.22744901424408D+04, 0.78779333020687D+03, &
    -69.830527374994D+00,  0.17863832875422D+05,-0.39514731563338D+05, &
    0.33803884280753D+05, -0.13855050202703D+05,-0.25637436613260D+06, &
    0.48212575981415D+06, -0.34183016969660D+06, 0.12223156417448D+06, &
    0.11797433655832D+07, -0.21734810110373D+07, 0.10829952168620D+07, &
   -0.25441998064049D+06, -0.31377774947767D+07, 0.52911910757704D+07, &
   -0.13802577177877D+07, -0.25109914369001D+06, 0.46561826115608D+07, &
   -0.72752773275387D+07,  0.41774246148294D+06, 0.14016358244614D+07, &
   -0.31555231392127D+07,  0.47929666384584D+07, 0.40912664781209D+06, &
   -0.13626369388386D+07,  0.69625220862664D+06,-0.10834900096447D+07, &
   -0.22722827401688D+06,  0.38365486000660D+06, 0.68833257944332D+04, &
    0.21757245522644D+05, -0.26627944829770D+04,-0.70730418082074D+05, &
   -0.225D+00, -1.68D+00, 0.055D+00, -93.0D+00 /)
  real gr
  real hr
  integer i
  integer, parameter, dimension ( 40 ) :: ii = (/ &
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6, &
    8,8,8,8,2,2,0,4,2,2,2,4 /)
  integer j
  integer, parameter, dimension ( 40 ) :: jj = (/ &
    2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,&
    2,3,5,7,1,4,4,4,0,2,0,0 /)
  integer k
  integer l
  integer nc
  real pr
  real q10
  real q20
  real q2a
  real q5t
  real qm
  real qp
  real qr(11)
  real qt(10)
  real rho
  real sr
  real t
  real tau
  real tx
  real ur
  real v
!
  ERRTOL = sqrt ( epsilon ( ERRTOL ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'RESID - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'RESID - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if

  nc = 36
  dpdrr = 0.0D+00
  pr = 0.0D+00
  ar = 0.0D+00
  dadt = 0.0D+00
  cvr = 0.0D+00
  dpdtr = 0.0D+00

  ex0 = - aa * rho
  ex0 = max ( ex0, -225.0 )
  ex0 = min ( ex0, 225.0 )
  e = exp ( ex0 )

  q10 = rho * rho * e
  q20 = 1.0D+00 - e

  qr(1) = 0.0D+00
  qr(2) = q10
  do i = 2, 10
    qr(i+1) = qr(i) * q20
  end do

  v = TZERO / t
  qt(1) = t / TZERO
  do i = 2, 10
    qt(i) = qt(i-1) * v
  end do

  do i = 1, nc

    k = ii(i) + 1
    l = jj(i)
    qp = g(i) * aa * qr(k+1) * qt(l+1)
    pr = pr + qp

    dpdrr = dpdrr + aa * ( 2.0D+00 / rho - ( 1.0D+00 - e * dble ( k - 1 ) / &
      ( 1.0D+00 - e ) ) ) * qp

    ar = ar + g(i) * qr(k+2) * qt(l+1) / ( rho**2 * e * dble ( k ) &
      * calc_gascon ( ) * t )

    dfdt = ( 1.0D+00 - e )**k * dble ( 1 - l ) * qt(l+2) / TZERO / dble ( k )

    dadt = dadt + g(i) * dfdt

    dpdtr = dpdtr + g(i) * dfdt * rho**2 * e * dble ( k ) / ( 1.0D+00 - e )

    cvr = cvr + g(i) * dble ( l ) * dfdt / calc_gascon()

  end do

  qp = 0.0D+00
  q2a = 0.0D+00
  do j = 37, 40
    k = ii(j)
    ddz = adz(j-36)
    del = rho / ddz - 1.0D+00

    if ( abs ( del ) < ERRTOL ) then
      del = ERRTOL
    end if

    ex1 = - aad(j-36) * del**k
    ex1 = max ( ex1, - 225.0 )
    ex1 = min ( ex1, 225.0 )
    dex = exp ( ex1 ) * del**jj(j)

    att = aat(j-36)
    tx = atz(j-36)
    tau = ( t / tx ) - 1.0D+00

    ex2 = - att * tau**2
    ex2 = max ( ex2, - 225.0 )
    ex2 = min ( ex2, 225.0 )
    q10 = dex * exp ( ex2 )

    qm = dble ( jj(j) ) / del - dble ( k ) * aad(j-36) * del**(k-1)
    fct = qm * rho**2 * q10 / ddz

    q5t = fct * ( 2.0D+00 / rho + qm / ddz ) - ( rho / ddz )**2 * q10 * &
      ( dble ( jj(j) ) / del**2 + dble ( k * ( k - 1 ) ) * aad(j-36) * &
      del**(k-2) )

    dpdrr = dpdrr + q5t * g(j)
    qp = qp + g(j) * fct
    dadt = dadt - 2.0D+00 * g(j) * att * tau * q10 / tx
    dpdtr = dpdtr - 2.0D+00 * g(j) * att * tau * fct / tx

    q2a = q2a + t * g(j) * att * ( 4.0D+00 * ex2 + 2.0D+00 ) * q10 / tx**2

    ar = ar + q10 * g(j) / ( calc_gascon() * t )

  end do

  cvr = cvr + q2a / calc_gascon()
  pr = pr + qp
  sr = - dadt / calc_gascon()
  ur = ar + sr
!
!  Assign dimensions.
!
  ar =  calc_gascon() * t *  ar
  cvr = calc_gascon() *     cvr
  sr =  calc_gascon() *      sr
  ur =  calc_gascon() * t *  ur
!
!  Adjust energies.
!
  ar = ar + calc_gascon ( ) * t * SREF - calc_gascon ( ) * UREF
  sr = sr - calc_gascon ( ) * SREF
  ur = ur - calc_gascon ( ) * UREF

  gr = ar + pr / rho
  hr = ur + pr / rho

  return
end subroutine resid

subroutine surten ( t, sigma )

!*****************************************************************************80
!
!! SURTEN returns the surface tension as a function of temperature.
!
!  Discussion:
!
!    SURTEN uses an equation that yields values of the surface tension to
!    within the accuracy of measurements from the triple point to the
!    critical point.
!
!      Sigma = B * ( (TSTAR-T)/TSTAR)**Mu * (1+b*(TSTAR-T)/TSTAR)
!
!    where:
!
!      TSTAR = 647.15 Degrees Kelvin,
!      B = 0.2358 Pascals * Meters
!      b = -0.625,
!      Mu = 1.256.
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) SIGMA, the surface tension,
!    in Pascal * m = Newton / m.
!
  implicit none

  real ( kind = 8 ), parameter :: b_cap = 0.2358D+00
  real ( kind = 8 ), parameter :: b_small = -0.625D+00
  real ( kind = 8 ), parameter :: mu = 1.256D+00
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_star = 647.15D+00
  real ( kind = 8 ) term
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SURTEN - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  term = ( t_star - t ) / t_star
  sigma = b_cap * term**mu * ( 1.0D+00 + b_small * term )
!
!  Need this conversion to match the table, but justification is there none.
!
  sigma = 1000.0D+00 * sigma

  return
end subroutine

subroutine tdpsdt ( t, dp )
!
!*******************************************************************************
!
!! TDPSDT computes the quantity T * dP(Sat)/dT.
!
!
!  Discussion:
!
!    Here T is the temperature and P(Sat) is the vapor pressure.  
!    It is used by TSAT_EST and TSAT.
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Output, real DP, the value T*(dP(Sat)/dT), 
!    in MegaPascals.
!
  real, parameter :: TREF = 647.25D+00
!
  real, parameter, dimension ( 8 ) :: a = (/ &
      -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,   -148.39348D+00,  48.631602D+00 /)
  real b
  real c
  real dp
  integer i
  real q
  real t
  real v
  real w
  real y
  real z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'TDPSDT - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if

  v = t / TREF
  w = 1.0D+00 - v
  b = 0.0D+00
  c = 0.0D+00
  do i = 1, 8
    z = dble ( i + 1 ) / 2.0D+00
    y = a(i) * w**z
    c = c + ( y / w ) * ( 0.5D+00 - 0.5D+00 * dble ( i ) - 1.0D+00 / v )
    b = b + y
  end do

  q = b / v
  dp = 22.093D+00 * exp ( q ) * c

  return
end subroutine tdpsdt

subroutine thercon ( t, rho, lambda )
!
!*******************************************************************************
!
!! THERCON calculates the thermal conductivity for given temperature and density.
!
!
!  Modified:
!
!    20 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the density, in G/CM3.
!
!    Output, real LAMBDA, the thermal conductivity, 
!    in mW/(m degrees Kelvin).
!
  real, parameter :: ACON = 18.66D+00
  real, parameter :: BCON = 1.00D+00
  real, parameter :: CCON = 3.7711D-08
  real, parameter :: OMEGA = 0.4678D+00
  real, parameter :: PREF = 22.115D+00
  real, parameter :: RHOREF = 317.763D+00
  real, parameter :: TREF = 647.27D+00
!
  real a
  real, parameter, dimension ( 0:3 ) :: acof = (/ &
    2.02223D+00, 14.11166D+00, 5.25597D+00, -2.01870D+00 /)
  real b(0:4,0:5)
  real chi
  real cjth
  real cjtt
  real cp
  real cv
  real dpdr
  real dpdr2
  real dpdt
  real eta
  real g
  real h
  integer i
  integer j
  real lambda
  real lambda0
  real lambda_del
  real p
  real power
  real rho
  real rho2
  real s
  real sum
  real t
  real u
!
  data b / &
    1.3293046D+00,    1.7018363D+00,   5.2246158D+00,  & 
    8.7127675D+00, -1.8525999D+00, &
   -0.40452437D+00,  -2.2156845D+00, -10.124111D+00,   &
   -9.5000611D+00,  0.93404690D+00, &
    0.24409490D+00,   1.6511057D+00,   4.9874687D+00, &  
    4.3786606D+00,  0.0D+00, &
    0.018660751D+00, -0.76736002D+00, -0.27297694D+00, &
    -0.91783782D+00, 0.0D+00, &
    -0.12961068D+00,  0.37283344D+00, -0.43083393D+00,  &
     0.0D+00,        0.0D+00, &
    0.044809953D+00, -0.11203160D+00,  0.13333849D+00,  &
    0.0D+00,        0.0D+00 /
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'THERCON - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'THERCON - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute DPDR, DPDT, ETA.
!
  call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

  call viscos ( t, rho, eta )
!
!  Convert RHO from G/CM3 to KG/M3,
!  Convert DPDR from ? to ?.
!
  rho2 = 1000.0D+00 * rho
  dpdr2 = dpdr / 1000.0D+00
!
!  Compute LAMBDA0.
!
  sum = 0.0D+00
  do i = 0, 3
    sum = sum + acof(i) * ( TREF / t )**i
  end do

  lambda0 = sqrt ( t / TREF ) / sum
!
!  Compute CHI.
!
  chi = rho2 * PREF / ( RHOREF**2 * dpdr2 )
!
!  Compute delta_Lambda
!
  power = - ACON * ( ( TREF - t ) / t )**2 - BCON * ( ( rho2 - RHOREF ) &
    / RHOREF )**4

  lambda_del = ( CCON / eta ) * ( ( t * RHOREF ) / ( TREF * rho ) )**2 &
    * ( TREF / PREF )**2 * dpdt**2 * chi**OMEGA * sqrt ( rho2 / RHOREF ) &
    * exp ( power )
!
!  Compute LAMBDA.
!
  sum = 0.0D+00
  do i = 0, 4
    do j = 0, 5
      sum = sum + b(i,j) * ( ( TREF - t ) / t )**i * &
        ( ( rho2 - RHOREF ) / RHOREF )**j
    end do
  end do

  lambda = lambda0 * exp ( ( rho2 / RHOREF ) * sum ) + lambda_del
!
!  Temporary fix.
!
  lambda = 1000.0D+00 * lambda

  return
end subroutine thercon

subroutine therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )
!
!*******************************************************************************
!
!! THERM calculates thermodynamic functions given temperature and density.
!
!
!  Discussion:
!
!    Thermodynamic values were calculated from an analytic equation
!    that approximates the Helmholtz function (specific Helmholtz
!    energy) for ordinary water and steam, of the form A=A(rho,T)
!    where A is the Helmholtz function, rho the density, and T
!    the absolute (thermodynamic) temperature.  Any thermodynamic
!    value for any state, liquid, vapor or metastable, may be
!    calculated by differentiation of this equation in accord with
!    the first and second laws of thermodynamics.
! 
!    The International Association for the Properties of Steam
!    has provisionally accepted this formulation for the range
!    273.15 <= T <= 1273.15 degrees Kelvin, where, for 423.15 <= T,
!    the maximum pressure is Pmax = 1500 MPa = 15000 bar, and for
!    273.15 <= T < 423.15, the maximum pressure is
!    Pmax = 100 * (5 + (T-273.15)/15) MPa.
! 
!    Close to the critical point, a small region is excluded:
!    Abs(T-Tk) < 1, abs((rho-rhok)/rhok) < 0.3.
! 
!    The equation has a wider useful range, namely, fluid states
!    of pure, undissociated water and steam defined by
!    260 <= T <= 2500 K and 0 <= P <= 3000 MPa.
! 
!    Thermodynamic property values for specific volume, density,
!    specific internal energy, specific enthalpy, and specific
!    entropy of water and steam were tabulated over the range
!    0 <= t <= 2000 C, 0 <= P <= 3000 MPa.  The reference
!    state is the liquid at the triple point, for which the
!    internal energy and entropy have been assigned the value zero.
!     
!    Thermodynamic quantities are determined from the Helmholtz function
!    A(rho,T), which is computed as the sum of three terms:
!
!      A(rho,T) = A(base)(rho,T) + A(residual)(rho,T) + A(ideal)(T)
!                                                       (Equation 1)
!
!    Because A(rho,T) is everywhere single valued and analytic,
!    we can derive closed form relations for all other properties.
!    In the following, unless otherwise indicated, the independent
!    variables are temperature T and density RHO, and differentiation
!    with respect to one variable is to imply that the other is fixed.
! 
!    Pressure:                  P       = RHO**2 * dA/dRHO
!    Density derivative:        dP/dRHO = 2*P/RHO + RHO**2 * d2A/dRHO2
!    Temperature derivative:    dP/dT   = RHO**2 * d2A/(dRHO dT)
!    Specific entropy:          S       = - dA/dT
!    Specific internal energy:  U       = A + T*S
!    Specific enthalpy:         H       = U + P/RHO
!    Specific heat capacity
!      at constant volume:      Cv      = - T * d2A/dT2
!    Specific Gibbs function:   G       = A + P/RHO
!    Specific heat capacity
!      at constant pressure:    Cp      = Cv + (T*(dP/dT)**2)/(RHO**2*dP/dRHO)
!    Speed of sound:            Omega   = Sqrt ((Cp/Cv) * dP/dRHO)
!    Second virial coefficient: B       = 1/(2*R*T) * (d2P/dRHO2) (at RHO=0)
!    Isothermal Joule-Thomson
!      coefficient:             DeltaT  = (dH/dP) (fixed T) =
!                                         (1/RHO)-(T*dP/dT)/(RHO**2*dP/dRHO)
!    Joule-Thomson coefficient: Mu      = (dT/dP) (fixed H) = DeltaT/Cp
!    Isentropic temperature-
!      pressure coefficient:    BetaS   = (dT/dP) (fixed S) =
!                                         (DeltaT - 1/RHO)/Cp
!   
!  Modified:
!
!    19 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the fluid density, in G/CM3.
!
!    Output, real A, the Helmholtz function, in KJ/KG.
!
!    Output, real CJTH, the Joule-Thomson coefficient,
!    in K/MegaPascals.
!
!    Output, real CJTT, the isothermal Joule-Thomson coefficient,
!    in CM3/G.
!
!    Output, real CP, the isobaric (constant pressure) heat
!    capacity, in KJ/(KG degrees Kelvin).
!
!    Output, real CV, the isochoric (constant volume) heat capacity,
!    in KJ/(KG degrees Kelvin).
!
!    Output, real DPDR, the partial derivative 
!    dP(T,RHO)/dRHO, with T held fixed, in MegaPascals*CM3/G.
!
!    Output, real DPDT, the partial derivative 
!    dP(T,RHO)/dT, with RHO held fixed, in MegaPascals/degrees Kelvin.
!
!    Output, real G, the Gibbs free energy, in KJ/KG.
!
!    Output, real H, the enthalpy, in KJ/KG.
!
!    Output, real P, the pressure, in MegaPascals.
!
!    Output, real S, the entropy, in KJ/(KG degrees Kelvin).
!
!    Output, real U, the internal energy, in KJ/KG.
!
  real a
  real ab
  real ai
  real ar
  real cjth
  real cjtt
  real cp
  real cpi
  real cv
  real cvb
  real cvi
  real cvr
  real dpdr
  real dpdrb
  real dpdrr
  real dpdt
  real dpdtb
  real dpdtr
  real g
  real gb
  real gi
  real gr
  real h
  real hb
  real hi
  real hr
  real p
  real pb
  real pr
  real rho
  real s
  real sb
  real si
  real sr
  real t
  real u
  real ub
  real ui
  real ur
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'THERM - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'THERM - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if

  call ideal ( t, ai, cpi, cvi, gi, hi, si, ui )

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  a =       ab +    ar +  ai
  cv =     cvb +   cvr + cvi
  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr
  p =       pb +    pr
  s =       sb +    sr +  si
  u =       ub +    ur +  ui

  g = a + p / rho
  h = u + p / rho
  cp = cv + t * dpdt**2 / ( dpdr * rho**2 )
  cjtt = 1.0D+00 / rho - t * dpdt / ( dpdr * rho**2 )
  cjth = - cjtt / cp

  return
end subroutine therm

subroutine tsat ( p, t, rhol, rhov )
!
!*******************************************************************************
!
!! TSAT calculates the saturation temperature for a given pressure.
!
!
!  Discussion:
!
!    The corresponding liquid and vapor densities are also computed.
!    The saturation temperature is also known as the "vapor temperature".
!
!  Modified:
!
!    21 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real P, the vapor pressure, in MegaPascals.
!
!    Output, real T, the vapor temperature, in degrees Kelvin.
!
!    Input/output, real RHOL, the liquid density, in G/CM3.
!    On input, if RHOL is nonnegative, it is used as an estimate for
!    the density.
!
!    Input/output, real RHOV, the vapor density, in G/CM3.
!    On input, if RHOV is nonnegative, it is used as an estimate for
!    the density.
!
  integer, parameter :: ITMAX = 50
!
  real delg
  real dp
  real dp_BOT
  real dp2
  real ERRTOL
  integer it
  real p
  real rhol
  real rhov
  real t
!
  ERRTOL = sqrt ( epsilon ( ERRTOL ) )
!
!  Refuse to handle zero or negative pressure.
!
  if ( p <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'TSAT - Fatal error!'
    write(*, * ) '  The input pressure must be positive!'
    write(*, * ) '  Your value was P = ', p
    stop
  end if
!
!  Estimate the saturation temperature.
!
  call tsat_est ( p, t )

  do it = 1, ITMAX

    call corr ( t, p, rhol, rhov, delg )

    dp_bot = 1.0D+00 / rhov - 1.0D+00 / rhol
    dp_bot = max ( dp_bot, ERRTOL )

    dp = delg * calc_gascon ( ) * t / dp_bot

    call tdpsdt ( t, dp2 )

    t = t * ( 1.0D+00 - dp / dp2 )

    if ( abs ( delg ) < ERRTOL ) then
      return
    end if

  end do

  write(*, * ) ' '
  write(*, * ) 'TSAT - Warning!'
  write(*, * ) '  The iteration did not converge.'
  write(*, * ) '  Number of iterations was ', ITMAX
  write(*, * ) '  Last iterate was ', t
  write(*, * ) '  Last DELG was ', delg
   
  return
end subroutine tsat

subroutine tsat_est ( p, t )
!
!*******************************************************************************
!
!! TSAT_EST makes a rough estimate of the saturation temperature.
!
!
!  Discussion:
!
!    The saturation temperature is also called the vapor temperature.
!
!  Modified:
!
!    22 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real P, the pressure, in MegaPascals.  The tabulated
!    range for P is 0.00061173 MegaPascals <= P <= 22.055 MegaPascals.  
!    The input value of P must be positive.
!
!    Output, real T, the saturation temperature, 
!    in degrees Kelvin.
!
  integer, parameter :: KMAX = 8
  real, parameter :: PMAX = 22.055D+00
  real, parameter :: TMAX = 647.126D+00
  real, parameter :: TMIN = 273.15D+00
!
  real dp
  real dt
  real ERRTOL
  integer k
  real p
  real pl
  real pp
  real t
  real told
!
  ERRTOL = sqrt ( epsilon ( ERRTOL ) )
!
!  Refuse to handle zero or negative pressure.
!
  if ( p <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'TSAT_EST - Fatal error!'
    write(*, * ) '  The input pressure must be positive!'
    write(*, * ) '  Your value was P = ', p
    stop
  end if

  if ( p > PMAX ) then
    t = TMAX
    return
  end if
 
  pl = 2.302585D+00 + log ( p )

  t = 372.83D+00 &
    + pl * ( 27.7589D+00 &
    + pl * ( 2.3819D+00 &
    + pl * ( 0.24834D+00 &
    + pl *   0.0193855D+00 ) ) )

  t = min ( t, TMAX )
  t = max ( t, TMIN )

  dt = 0.0D+00

  do k = 1, KMAX

    call psat_est ( t, pp )

    call tdpsdt ( t, dp )

    if ( abs ( p - pp ) < ERRTOL * p ) then
      return
    end if

    dt = t * ( p - pp ) / dp
 
    told = t
    t = t * ( 1.0D+00 + ( p - pp ) / dp )
    t = min ( t, TMAX )
    t = max ( t, TMIN )

    if ( abs ( dt ) < ERRTOL * ( abs ( t ) + 1.0D+00 ) ) then
      return
    else if ( abs ( t - told ) < ERRTOL ) then
      return
    end if

  end do

  write(*, * ) ' '
  write(*, * ) 'TSAT_EST - Warning!'
  write(*, * ) '  The iteration did not converge.'
  write(*, * ) '  Number of iterations was ', KMAX
  write(*, * ) '  Convergence tolerance was ', ERRTOL
  write(*, * ) '  Last iterate was ', t
  write(*, * ) '  Last correction was ', dt

  return
end subroutine tsat_est

subroutine viscos ( t, rho, eta )
!
!*******************************************************************************
!
!! VISCOS calculates the viscosity for given temperature and density.
!
!
!  Verification:
!
!    This subroutine does not produce the values in the tables.
!    It is very badly off.
!
!    In a discussion by John Burkardt posted on his personal website...
!
!    "On 02 February 2002, I discovered that the Haar/Gallagher/Kell
!    reference apparently reversed the sign on the A3 coefficient.
!    That made the results better, but still off.
!
!    Apparently Haar/Gallagher/Kell had a transcription error in
!    the value of B(4,1), which they list as -0.273093, but which
!    should be -0.253093.  
!
!    These two corrections courtesy of Meyer/McClintock/Silvestri/Spencer.
!
!    Now the results look proper!  And just 12 years late..."
!
!  Modified:
!
!    20 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    International Association for the Properties of Steam,
!    Release on Dynamic Viscosity of Water Substance,
!    National Bureau of Standards, Washington DC, 1975, revised 1983.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the density, in G/CM3.
!
!    Output, real ETA, the viscosity, in MegaPascal seconds.
!
  real, parameter :: RHOREF = 0.317763D+00
  real, parameter :: TREF = 647.27D+00
!
  real, parameter, dimension ( 0:3 ) :: a = (/ &
    0.0181583D+00, 0.0177624D+00, 0.0105287D+00, -0.0036744D+00 /)
  real b(0:5,0:4)
  logical, parameter :: LOCAL_DEBUG = .false.
  real eta
  real eta0
  integer i
  integer j
  integer k
  real rho
  real sum
  real t
!
  data b / &
    0.501938D+00,   0.162888D+00,  -0.130356D+00,  &
    0.907919D+00, -0.551119D+00,   0.146543D+00,   &
    0.235622D+00,   0.789393D+00,   0.673665D+00,  &
    1.207552D+00,  0.0670665D+00, -0.0843370D+00,  &
   -0.274637D+00,  -0.743539D+00,  -0.959456D+00,  &
   -0.687343D+00, -0.497089D+00,   0.195286D+00,   &
    0.145831D+00,   0.263129D+00,   0.347247D+00,  &
    0.213486D+00,  0.100754D+00,  -0.032932D+00,   &
   -0.0270448D+00, -0.0253093D+00, -0.0267758D+00, &
   -0.0822904D+00, 0.0602253D+00, -0.0202595D+00 /
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'VISCOS - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'VISCOS - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute ETA0.
!
  sum = 0.0D+00
  do k = 0, 3
    sum = sum + a(k) * ( TREF / t )**k
  end do

  eta0 = sqrt ( t / TREF ) / sum

  if ( LOCAL_DEBUG ) then
    write(*, * ) 'VISCOS : ETA0 = ', eta0
  end if
!
!  Compute ETA.
!
  sum = 0.0D+00
  do i = 0, 5
    do j = 0, 4
      sum = sum + b(i,j) * ( ( TREF - t ) / t )**i * &
        ( ( rho - RHOREF ) / RHOREF )**j
    end do
  end do

  eta = eta0 * exp ( ( rho / RHOREF ) * sum )

  return
end subroutine viscos

subroutine calc_volume ( t, rho, v, dvdt, dvdr )
!
!*******************************************************************************
!
!! VOLUME computes specific volume derivatives given temperature and density.
!
!
!  Discussion:
!
!    Because A(rho,T) is everywhere single valued and analytic,
!    we can derive closed form relations for all other properties.
!    In the following, unless otherwise indicated, the independent
!    variables are temperature T and density RHO, and differentiation
!    with respect to one variable is to imply that the other is fixed.
! 
!  Modified:
!
!    28 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real T, the temperature, in degrees Kelvin.
!
!    Input, real RHO, the fluid density, in G/CM3.
!
!    Output, real V, the specific volume, in CM3/G.
!
!    Output, real DVDT, the partial derivative dV(T,RHO)/dT, 
!    where V is the specific volume, in CM3 / (G * degrees Kelvin).
!
!    Output, real DVDT, the partial derivative dV(T,RHO)/dRHO, 
!    where V is the specific volume, in CM3**2 / ( G**2 ).
!
  real ab
  real ar
  real cvb
  real cvr
  real dpdr
  real dpdrb
  real dpdrr
  real dpdt
  real dpdtb
  real dpdtr
  real dvdr
  real dvdt
  real gb
  real gr
  real hb
  real hr
  real pb
  real pr
  real rho
  real sb
  real sr
  real t
  real ub
  real ur
  real v
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'VOLUME - Fatal error!'
    write(*, * ) '  Input temperature T must be positive.'
    write(*, * ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write(*, * ) ' '
    write(*, * ) 'VOLUME - Fatal error!'
    write(*, * ) '  Input density RHO must be positive.'
    write(*, * ) '  Input value was RHO = ', rho
    stop
  end if

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr

  dvdt = dpdt / ( dpdr * rho**2 )
  dvdr = - 1.0D+00 / rho**2
  v = 1.0D+00 / rho

  return
end subroutine calc_volume
!

subroutine calchfg(hfg,p)
!  This subroutine calculates hfg (a.k.a. the latent heat of vapourisation) from the pressure
!  It uses the NBC/NRC Steam Table routines 
   implicit none
   real, intent(in) :: p     ! in MPa
   real, intent(out):: hfg   ! in kJ/kg or J/g
!  local variables
   real             :: ts,rhof,rhog,hf,hg
!  local variables - these are not actually used, they are just required for the call to therm
   real             :: a,cjth,cjtt,cp,cv,dpdr,dpdt,g,psat,s,u 
   
!  Determine saturation temperature, at the total pressure, P.  rhof and rhog are set as negative
!  so that they are NOT used as estimates of density when calculating Ts.  The tsat subroutine then
!  returns tsat, and the saturation densities rhof and rhog.
   rhof = -1.0
   rhog = -1.0
   call tsat (p, ts, rhof, rhog)

!  Now we have the saturation temperature and the saturation densities we can call therm.  Since 
!  combinations of temperature and density and thermodynamically single valued everywhere
!  then therm can calculate all of the other thermodynamic properties.  Here we call it to 
!  calculate the saturation enthalpies hf and hg. 
   
   call therm ( ts, rhof, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, hf, psat, s, u )
   call therm ( ts, rhog, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, hg, psat, s, u )
   
   hfg = hg - hf
   
end subroutine calchfg


!============================================================
!
SUBROUTINE SILICAEOS(PCOEF,SCOEF,DENST,PBF,PATMOS,TBF)
!
! This subroutine calculates the solid-liquid-vapor equilibrium
! of silica over a wide range of temperature and pressure. Indeed
! the following parameters should be remembered:
!   Tc = 11976.0 K ; Pc = 2.0 Kbar ; DENSc = 0.58 g/cm**3
!   Tmelt (P = 1 bar)        =   1996.0 K
!   Tglass (P = 1 bar)       =   1480.0 K
!   Tboil (P = 1 bar         =   2530.0 K
!   Tboil (Psat = 70 bar)    =   5638.0 K
!   DENST (g, T = 9000 K     =   0.019 g/cm**3
!   DENST (l, T = 1673.0 K)  =   2.20  g/cm**3
!   DENST (l, T = 4000.0 K)  =   2.50  g/cm**3
!   Sound Speed (P,T = room) =   5760.0 m/s
! We still have a lot of work on the blending of EOS used for the 
! SLVE of silica.  For the solid-liquid glass transition, a 
! Murnagham EOS is used as described by Mao et al. in Journal of 
! Alloys and Compounds 327(2001):253-262. For the liquid phase,
! the fitted equation suggested by Guissani & Guillot, J. Chem.
! Phys. 104(1996):7633-7644 is used.  Unfortunately, no reliable 
! EOS for the gas phase was found, therefore, a constant density
! was used.
!
  IMPLICIT NONE
  REAL::PBF,PATMOS,TBF,PCOEF,SCOEF
  REAL::T,P,TBOIL
  REAL,PARAMETER::TC=11976.0,RHOC=0.58,PC=2069.0
  REAL,PARAMETER::TMELT=1996.0,MW=60.08,DENSS=2.66,&
       SOUND=5760.0E+2
  REAL,PARAMETER::DELTA=50.0 ! Delta is applied to the coupled 
  ! temperature and density changings
  REAL,PARAMETER::TZERO=273.15,PCONV=1.0E-6,CONVSO=1.0E+2
  ! Converter coefficientes:   for C --> K and
  !                            for g/(cm.s**2) --> bar
  !                            for m/s --> cm/s
  REAL::RHOL,RHOL2,RHOG,DENST
  REAL::TAU,TAU2
  REAL::PREF,DENSTREF,TBOILREF,RHOGREF,RHOLREF
  !
! these are part of module now, so shouldn't declare them (same goes for
! EXPREP throughout this file...)
!  REAL::SIGOI,CALCBOIL,FCEOS1,FCEOS2,CALCSOUND,DSTPT ! Real functions used in the program

! Melting point is set as a constant, 1996 K, while the boiling
! point starts at 2530.15 K at 1 bar and is increased linearly
! up to 5638 K at 70 bar
! PARA APLICAR A PERTUBACAO PARA ENCONTRAR PCOEF, USE UMA FC PARA CALCULAR DENSIDADE
! A PARTIR DA DENSIDADE NORMAL COM INCREMENTOS EM PRESSAO E TEMPERATURA !!
!
      p=PBF+PATMOS
      PREF=PATMOS
      t=TBF
!
    TBOIL=CALCBOIL(P)
    TBOILREF=CALCBOIL(PREF) ! Reference state
!
     IF(T<TMELT-DELTA)THEN
        CALL BEFORMELT(P,T,DENST)
        CALL BEFORMELT(PREF,T,DENSTREF)
     ELSEIF((T>=TMELT-DELTA).AND.(T<=TMELT+DELTA))THEN
        CALL BEFORMELT(P,T,RHOG)
        CALL BEFORMELT(PREF,T,RHOGREF)
        TAU=ABS(1.0-T/TC)
        RHOL=(FCEOS1(TAU)+FCEOS2(TAU)+1.0)*RHOC 
        RHOLREF=DSTPT(RHOL,P,T,PREF,T)
        DENST=SIGOI(T,TMELT,DELTA,RHOL,RHOG)
        DENSTREF=SIGOI(T,TMELT,DELTA,RHOLREF,RHOGREF)
     ELSEIF((T>TMELT+DELTA).AND.(T<TBOIL-DELTA))THEN
        TAU=ABS(1.0-T/TC)
        RHOL=(FCEOS1(TAU)+FCEOS2(TAU)+1.0)*RHOC 
        RHOL=DSTPT(RHOL,P,T,(PBF+PATMOS)*PCONV,TBF+TZERO)
        RHOLREF=DSTPT(RHOL,P,T,PREF,T)
        RHOG=0.0
        DENST=RHOL
        DENSTREF=RHOLREF
     ELSEIF((T>=TBOIL-DELTA).AND.(T<=TBOIL+DELTA))THEN
        TAU=ABS(1.0-T/TC)
        RHOL=(FCEOS1(TAU)+FCEOS2(TAU)+1.0)*RHOC 
        RHOL=DSTPT(RHOL,P,T,(PBF+PATMOS)*PCONV,TBF+TZERO)
        RHOLREF=DSTPT(RHOL,P,T,PREF,T)
        DENST=SIGOI(T,TBOIL,DELTA,0.019,RHOL)
        DENSTREF=SIGOI(T,TBOILREF,DELTA,0.019,RHOLREF)
     ELSEIF(T>=TBOIL+DELTA)THEN!.AND.(T<8000.0))THEN
        DENST=0.019  ! This and the next 'elseif' may be changed later, when a 
                     ! reliable EOS for the gas phase is found.
        DENST=DSTPT(DENST,P,T,(PBF+PATMOS)*PCONV,TBF+TZERO)
        DENSTREF=DSTPT(DENST,P,T,PREF,T)
     ENDIF
     PCOEF=ABS((DENST-DENSTREF)/MAX(1.0E-5,ABS(P/PCONV-PREF/PCONV)))
     SCOEF=DENST-PCOEF*(PBF+PATMOS)*POWERFN(PCONV,-1.0)
 
!
  RETURN 
!
END SUBROUTINE SILICAEOS


!===================================================================
REAL FUNCTION CALCBOIL(P)
! This function calculates the boiling temperature from a 
! given pressure.  As I did not find any simple relationship
! between them (except the proper minimization function which
! would be very time consuming) I am assuming a linar function
! as: TBOIL= CONV*(P-1)+TBOIL0, where TBOIL0 is the boiling 
! temperature at 1 bar. However, I will assume that beyond 
! 70 bars, the boiling temperature remains constant.
!
!  use FLDebug
  IMPLICIT NONE
  REAL::P
  REAL,PARAMETER::CONV=45.0451304,TBOIL0=2530.15
!
    IF(P<=70.0)THEN
      CALCBOIL=(P-1.0)*CONV+TBOIL0
    ELSE
      CALCBOIL=69.0*CONV+TBOIL0
    ENDIF
!
  RETURN
!
END FUNCTION CALCBOIL

!=================================================================
!
!
REAL FUNCTION FCEOS1(TAU)
!
  IMPLICIT NONE
  REAL::TAU,EXPLOG
!
    EXPLOG=1.0/(1.0+EXP(-6.0*(1.2+LOG(TAU))))
    FCEOS1=4.05*POWERFN(TAU,0.5)*(0.62*POWERFN(TAU,-0.182)-0.25*&
           POWERFN(TAU,7.0))*EXPLOG 
!
  RETURN
END FUNCTION FCEOS1

!================================================================
!
REAL FUNCTION FCEOS2(TAU)
  IMPLICIT NONE
  REAL::TAU,EXPLOG
!
    EXPLOG=1.0/(1.0+EXP(-6.0*(1.8*LOG(TAU))))
    FCEOS2=1.0+7.50*TAU*(0.20*POWERFN(TAU,-0.45)-0.13*&
           POWERFN(TAU,7.0))*EXPLOG 
!
  RETURN
END FUNCTION FCEOS2

!================================================================
!
SUBROUTINE BEFORMELT(P,T,DENST)
  IMPLICIT NONE
  REAL::T,P,V,DENST
!  REAL,PARAMETER,DIMENSION(2)::BETAS=(/ 24.8361E-12,&
!          379.561E-17 /),BETAL=(/ 1.15032E-12,&
!          -114.60E-17 /)
  REAL,PARAMETER,DIMENSION(2)::BETAS=(/ 1.15032E-12,&
          -114.6E-17 /),BETAL=(/ 1.15032E-12,&
          -114.60E-17 /)
  REAL::KT,KP
  REAL,PARAMETER::KPS=6.0,KPL=7.923
!  REAL::V1T  - part of module, so don't declare
  LOGICAL::OPT
!
    OPT=.TRUE.
    IF(OPT)THEN
      KT=1.0/(BETAS(1)+BETAS(2)*T)
      KP=KPS
    ELSE    
      KT=1.0/(BETAL(1)+BETAL(2)*T)
      KP=KPL
    ENDIF
!
    V=V1T(T,OPT)*POWERFN(1.0+KP*P/KT,-1.0/KP)
!
    DENST=1.0/(V/60.08)
!
  RETURN
!
END SUBROUTINE BEFORMELT

!==================================================================
!
REAL FUNCTION V1T(T,OPT)
!
  IMPLICIT NONE
  REAL::T
  LOGICAL::OPT
  REAL,PARAMETER,DIMENSION(2)::ALPHAS=(/ -1.15006E-5,&
          99.5453E-11 /),ALPHAL=(/ 19.0134E-5,&
          -10162.7E-11 /)
!  REAL,PARAMETER,DIMENSION(2)::ALPHAS=(/ 5.96296E-5,&
!          6.72697E-11 /),ALPHAL=(/ 19.0134E-5,&
!          -10162.7E-11 /)
  REAL,PARAMETER::TR=298.15
  REAL::ALPHA,V1TR,DIFFT
!
    DIFFT=T-TR
    IF(OPT)THEN
      ALPHA=ALPHAS(1)*(DIFFT)+0.5*ALPHAS(2)*POWERFN(DIFFT,2.0)
      V1TR=22.6357
    ELSE
      ALPHA=ALPHAL(1)*(DIFFT)+0.5*ALPHAL(2)*POWERFN(DIFFT,2.0)
      V1TR=23.1493
    ENDIF
!
    IF(ALPHA>1000.0)ALPHA=1000.0
    V1T=V1TR*EXP(ALPHA)
!
  RETURN
!
END FUNCTION V1T


!

!============================================================
!
REAL FUNCTION POWERFN(A,B)
!
! This fucntion calculates the power through the LOG.
!          K=A**B ==> K=EXP(B*LN(A))
!
  IMPLICIT NONE
  REAL::A,B
!  REAL::EXP
  !
    IF(A==0.0)THEN
       POWERFN=0.0
    ELSE
       POWERFN=EXP(B*LOG(ABS(A)))
    ENDIF
  !
    IF((A<=0.0).AND.(MOD(B,2.0)/=0))POWERFN=-POWERFN
  !
  RETURN
END FUNCTION POWERFN



!=================================================================
!
REAL FUNCTION SIGOI(Y,Y0,WIDTH,LOWMAG,UPRMAG)
!
! The sigmoid function, varies between (LOWMAG,UPRMAG).
! Y is the variable of the function, Y0 is the centre
! of the function and WIDTH is the width of the sigmoid
! function.  The function looks like:
!
!       --------
!               /
!                /
!                --------
!
  IMPLICIT NONE
  REAL::Y,Y0,WIDTH,LOWMAG,UPRMAG
  REAL::ALPHA!,EXPREP
!
    IF(Y-Y0<-0.9*WIDTH)THEN
      SIGOI=UPRMAG
    ELSEIF(Y-Y0>0.9*WIDTH)THEN
      SIGOI=LOWMAG
    ELSE
      ALPHA=20.0/(WIDTH)
      SIGOI=(UPRMAG-LOWMAG)/(1.0+EXP(ALPHA*(Y-Y0)))+LOWMAG
    ENDIF
!
  RETURN
!
END FUNCTION SIGOI

!==================================================================
!
REAL FUNCTION DSTPT(DENST,P,T,P0,T0)
  !
  IMPLICIT NONE
  REAL::DENST,P,T,P0,T0
  REAL::AUX1,AUX2!EXPREP,
  !      
  IF(ABS(P-P0)<=1.0E-5)THEN
     AUX1=0.0
  ELSE
     AUX1=POWERFN(ABS(P-P0),-0.1)
     AUX1=EXP(-10.0*AUX1)
  ENDIF
  IF(ABS(T-T0)<=1.0E-5)THEN
     AUX2=0.0
  ELSE
     AUX2=POWERFN(ABS(T-T0),-1.0)
     AUX2=EXP(-273.15*AUX2)
  ENDIF
  DSTPT=DENST*(1.0+0.1*AUX1-0.01*AUX2)
RETURN
END FUNCTION DSTPT




!==================================================================
!
!==================================================================
!
!==================================================================
!
!  
!  
!
       SUBROUTINE NEUEQNSTA2(PCOEF,SCOEF,DENSTY,P,PATMOS,TEMP,EQNSTA)

       IMPLICIT NONE
       REAL PCOEF,SCOEF,TEMP,P,PATMOS,DENSTY
       INTEGER EQNSTA
! This sub finds the density DENSTY from TEMP and PRES
! TEMP is in Kelvin and density is in g/cm3
! NB DENSTY=PCOEF*(P+PATMOS)+SCOEF
! EQNSTA=7 lime EoS
! EQNSTA=8 Silica EoS
! EQNSTA=9 granite EoS
       REAL PRES,PRES2,DENST2
       PRES=P+PATMOS
       PRES2=1.0001*PRES
!
       CALL NEUEQNSTA(TEMP,PRES, DENSTY,EQNSTA)
       CALL NEUEQNSTA(TEMP,PRES2,DENST2,EQNSTA)
!  
       PCOEF=(DENSTY-DENST2)/(PRES-PRES2)
       SCOEF=DENSTY-PCOEF*PRES
       RETURN
       END SUBROUTINE NEUEQNSTA2
!
!
!
!
       SUBROUTINE NEUEQNSTA(TEMP,PRES,DENSTY,EQNSTA)

       IMPLICIT NONE
       REAL TEMP,PRES,DENSTY
       INTEGER EQNSTA
! This sub finds the density DENSTY from TEMP and PRES
!       INTEGER :: NLAYER1,NLAFIN,NTRAIN,NONODS,MXNCOLM,NLAYERS,NITS
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! EQNSTA=7 lime EoS
! EQNSTA=8 Silica EoS
! EQNSTA=9 granite EoS
       INTEGER, PARAMETER :: NLAYER1=2,NLAFIN=1,NONODS=15
!       INTEGER, PARAMETER :: NLAYER1=2,NLAFIN=1,NONODS=17,NTRAIN=10000
       INTEGER, PARAMETER :: MXNCOLM=NONODS*NONODS,NLAYERS=4
       REAL :: WEIGHT(MXNCOLM),WEINEW(MXNCOLM)
       REAL :: MINT,MAXT,MINP,MAXP,MIND,MAXD
       REAL :: FUNC,FUNNEW,PERTERB
       INTEGER :: COLM(MXNCOLM),FINDRM(NONODS+1)
       INTEGER :: NLAY(NLAYERS)
       INTEGER :: ITNGOOD,ITS,COUNT,ILAYER,I,NCOLM, IRED,ITRAIN,I1
       INTEGER :: LOPT,NOREXP
!
       NLAY(1)=2
       NLAY(2)=6
       NLAY(3)=6
       NLAY(4)=1
!
       LOPT=EQNSTA-6
!
! Define the sparcity of the matrix and set WEIGHT
       CALL SPAWEI(WEIGHT,MXNCOLM,NONODS,  &
                   NCOLM,FINDRM,COLM,NLAY,NLAYERS)
!
       CALL DEFNWEI(WEIGHT,NCOLM,LOPT,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
!
       CALL EOS_NEU(TEMP,PRES,DENSTY, &
           MINT,MAXT,MINP,MAXP,MIND,MAXD, &
           WEIGHT,NONODS, MXNCOLM, &
           NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
!
       RETURN
       END SUBROUTINE NEUEQNSTA
!
!      
!
!
       SUBROUTINE DEFNWEI(WEIGHT,NCOLM,LOPT,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,LOPT,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
       IF(LOPT.EQ.1) THEN
! Eqn of state for lime...
       CALL DEFNWEILIM(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ELSE IF(LOPT.EQ.2) THEN
! Eqn of state for Silica...
       CALL DEFNWEISIL(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ELSE IF(LOPT.EQ.3) THEN
! Eqn of state for Granite...
       CALL DEFNWEIGRA(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ENDIF
!
        RETURN
        END SUBROUTINE DEFNWEI
!
!
!
!
       SUBROUTINE DEFNWEILIM(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for lime...
        NOREXP=0 
!
        WEIGHT(1)=      5.97685132589016     
        WEIGHT(2)=      1.10400778174330     
        WEIGHT(3)=      -21.2825376649738     
        WEIGHT(4)=      0.169090329919894     
        WEIGHT(5)=       20.1788017885642     
        WEIGHT(6)=      0.131574759151893     
        WEIGHT(7)=       12.3755801537776     
        WEIGHT(8)=      0.766587143843660     
        WEIGHT(9)=      -24.4043122338779     
        WEIGHT(10)=      0.232113849596108     
        WEIGHT(11)=       4.85706354480158     
        WEIGHT(12)=      -1.26789757374056     
        WEIGHT(13)=     -0.255388642114362     
        WEIGHT(14)=       12.0554633479762     
        WEIGHT(15)=      -6.94915549914577     
        WEIGHT(16)=      -2.09264002488505     
        WEIGHT(17)=       15.8594342730431     
        WEIGHT(18)=       1.72679479920685     
        WEIGHT(19)=      -2.06771382299692     
        WEIGHT(20)=      0.395823993317512     
        WEIGHT(21)=      -1.12455553498827     
        WEIGHT(22)=      0.750060594828952     
        WEIGHT(23)=     -0.473464061154712     
        WEIGHT(24)=      -1.09279477569545     
        WEIGHT(25)=       4.10419973811868     
        WEIGHT(26)=      -15.2397251897475     
        WEIGHT(27)=       14.8805451851581     
        WEIGHT(28)=       9.99889550046277     
        WEIGHT(29)=      -20.6890317411573     
        WEIGHT(30)=       5.36436208014150     
        WEIGHT(31)=      -1.24892426108564     
        WEIGHT(32)=       1.51633332142468     
        WEIGHT(33)=      -3.14208255247642     
        WEIGHT(34)=      -2.56368066528082     
        WEIGHT(35)=       9.65149011111481D-002
        WEIGHT(36)=     -0.507349715070754     
        WEIGHT(37)=      0.423474661102974     
        WEIGHT(38)=      -29.6165259232046     
        WEIGHT(39)=       12.3201603741371     
        WEIGHT(40)=       4.00375892664630     
        WEIGHT(41)=      -33.2267781535008     
        WEIGHT(42)=      -3.26392133049933     
        WEIGHT(43)=      -1.24302679976193     
        WEIGHT(44)=       4.59939718630899     
        WEIGHT(45)=      -5.25970543903552     
        WEIGHT(46)=      -3.00135185108571     
        WEIGHT(47)=       6.34139347514880     
        WEIGHT(48)=      -1.84842529815566     
        WEIGHT(49)=      -18.5019808067015     
        WEIGHT(50)=      1.88379314700385     
        WEIGHT(51)=      17.5824625711418     
        WEIGHT(52)=    -0.527333559763088     
        WEIGHT(53)=     -18.5389144969952     
        WEIGHT(54)=      6.24702737659637 
!    MINT,MAXT,MINP,MAXP,MIND,MAXD...
        MINT=298.150000000000        
        MAXT=37340.0952636910        
        MINP=1013250.00000000     
        MAXP=126898713.821683        
        MIND=1.31000000000000        
        MAXD=3.34412549863268 
        RETURN
        END SUBROUTINE DEFNWEILIM
!
!
!
!
       SUBROUTINE DEFNWEISIL(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for Silica...
        NOREXP=0 
!
        WEIGHT(1)=   6.12626077572350     
        WEIGHT(2)=    0.272564908604740     
        WEIGHT(3)=    -9.75616050871453     
        WEIGHT(4)=    0.246912396109541     
        WEIGHT(5)=    -6.53574055759190     
        WEIGHT(6)=   -0.384054877796388     
        WEIGHT(7)=     -8.44375031552827     
        WEIGHT(8)=    -3.329123698799069E-002
        WEIGHT(9)=     -8.27742549631889     
        WEIGHT(10)=    -9.890057731467744E-002
        WEIGHT(11)=      9.73715337228838     
        WEIGHT(12)=    -0.240771751819873     
        WEIGHT(13)=     -4.49184124614962     
        WEIGHT(14)=      5.95516268582194     
        WEIGHT(15)=      2.82549172297567     
        WEIGHT(16)=    -1.894941028182234E-002
        WEIGHT(17)=     0.601560527201315     
        WEIGHT(18)=     -4.47260922589827     
        WEIGHT(19)=      2.09913243736491     
        WEIGHT(20)=     -20.5263926531534     
        WEIGHT(21)=     -12.2964085736977     
        WEIGHT(22)=     -15.5042442607703     
        WEIGHT(23)=     -17.9007234974105     
        WEIGHT(24)=      14.4442370983336     
        WEIGHT(25)=      1.31102712743980     
        WEIGHT(26)=      2.24711186705534     
        WEIGHT(27)=     0.862146059851389     
        WEIGHT(28)=      7.08865488987401     
        WEIGHT(29)=      4.18955246476026     
        WEIGHT(30)=    -0.989849174284737     
        WEIGHT(31)=     1.447487508202793E-002
        WEIGHT(32)=      3.48014767000874     
        WEIGHT(33)=      2.13069663254350     
        WEIGHT(34)=     -4.73426230781354     
        WEIGHT(35)=      4.49227690172869     
        WEIGHT(36)=      3.47395096346542     
        WEIGHT(37)=     -1.19391789986029     
        WEIGHT(38)=     -4.64835965448362     
        WEIGHT(39)=    -0.635582999915386     
        WEIGHT(40)=      2.15717035380372     
        WEIGHT(41)=     0.489452839441452     
        WEIGHT(42)=     -2.36516255020337     
        WEIGHT(43)=    -0.765529516551553     
        WEIGHT(44)=     -4.07263693916288     
        WEIGHT(45)=      2.99555734116876     
        WEIGHT(46)=     -3.86497423995451     
        WEIGHT(47)=     -2.35841013369974     
        WEIGHT(48)=     -1.43724880804879     
        WEIGHT(49)=      3.73463484818312     
        WEIGHT(50)=    -21.2831183705293     
        WEIGHT(51)=    0.384125180599244     
        WEIGHT(52)=     2.24339682643455     
        WEIGHT(53)=    -2.24290011447767     
        WEIGHT(54)=    0.146262916094321  
!    MINT,MAXT,MINP,MAXP,MIND,MAXD...   
   MINT=298.150000000000        
   MAXT=37340.0952636910        
   MINP=2026500.00000000     
   MAXP=127911963.821683       
   MIND=1.911301743049916E-002   
   MAXD=2.70021028300145 
   RETURN
   END SUBROUTINE DEFNWEISIL
!
!    
!
!
       SUBROUTINE DEFNWEIGRA(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for Granite...
        NOREXP=1
!
        WEIGHT(1)=     -1.25210991901872     
        WEIGHT(2)=      1.84286851967128D-002
        WEIGHT(3)=     -1.14594954358655     
        WEIGHT(4)=      8.40714881584379D-002
        WEIGHT(5)=      2.65597303171804D-002
        WEIGHT(6)=      1.55508856931543D-002
        WEIGHT(7)=     0.836002693041384     
        WEIGHT(8)=     -9.05909371896358D-002
        WEIGHT(9)=    -0.268794055455616     
        WEIGHT(10)=     -7.17034776815134D-002
        WEIGHT(11)=      1.24561087811569     
        WEIGHT(12)=     0.115202509920321     
        WEIGHT(13)=     0.592933066624559     
        WEIGHT(14)=     0.295443685199725     
        WEIGHT(15)=     0.170527010853434     
        WEIGHT(16)=     -1.60166187085996D-002
        WEIGHT(17)=     0.116064480970310     
        WEIGHT(18)=    -0.256250630064082     
        WEIGHT(19)=      1.67727052032752     
        WEIGHT(20)=      1.19632898578538     
        WEIGHT(21)=      1.82714810873958D-002
        WEIGHT(22)=    -0.791636710562007     
        WEIGHT(23)=     0.656387670046782     
        WEIGHT(24)=    -0.925300133084808     
        WEIGHT(25)=    -0.356175762842861     
        WEIGHT(26)=    -0.550524110962181     
        WEIGHT(27)=    -0.298790079514946     
        WEIGHT(28)=     0.322634049945166     
        WEIGHT(29)=     -2.40837140811772D-002
        WEIGHT(30)=      1.03851390805661     
        WEIGHT(31)=     0.616274464012693     
        WEIGHT(32)=     0.443745520350746     
        WEIGHT(33)=    -0.190649055120316     
        WEIGHT(34)=    -0.388445592379031     
        WEIGHT(35)=      9.01767363969453D-003
        WEIGHT(36)=    -0.385043937917201     
        WEIGHT(37)=     0.545501292960052     
        WEIGHT(38)=     0.472350823673588     
        WEIGHT(39)=    -0.232523211547267     
        WEIGHT(40)=    -0.392064393082278     
        WEIGHT(41)=    -0.116452739687640     
        WEIGHT(42)=    -0.516573848801829     
        WEIGHT(43)=     0.225262483882030     
        WEIGHT(44)=    -0.101955071126241     
        WEIGHT(45)=      5.47774259524076D-002
        WEIGHT(46)=     -9.01064058715637D-002
        WEIGHT(47)=    -0.101044266061742     
        WEIGHT(48)=     -9.15587454508020D-002
        WEIGHT(49)=     -8.02903612304057D-002
        WEIGHT(50)=     2.16295456236704     
        WEIGHT(51)=    -1.68201258314670     
        WEIGHT(52)=    0.415871026552132     
        WEIGHT(53)=    0.591087476532526     
        WEIGHT(54)=   -0.257449667903680     
!    MINT,MAXT,MINP,MAXP,MIND,MAXD... 
   MINT=298.150000000000        
   MAXT=37340.0952636910        
   MINP=2026500.00000000     
   MAXP=127911963.821683        
!   MIND=2276.12204764802        
!   MAXD=2675.75066083717        
   MIND=2.27612204764802        
   MAXD=2.67575066083717     
   RETURN
   END SUBROUTINE DEFNWEIGRA
!
!
!
!   
       SUBROUTINE SPAWEI(WEIGHT,MXNCOLM,NONODS, &
                   NCOLM,FINDRM,COLM,NLAY,NLAYERS)
! Define the sparcity of the matrix and set WEIGHT

       IMPLICIT NONE
       INTEGER :: MXNCOLM,NONODS,NCOLM,NLAYERS
       REAL :: WEIGHT(MXNCOLM)
       !!INTEGER :: NLAY(NLAYERS),FINDRM(NONODS+1),COLM(NCOLM)
       INTEGER :: NLAY(NLAYERS),FINDRM(NONODS+1),COLM(MXNCOLM)
! Local variables
       INTEGER :: NOD,COUNT,ILAYER,II,JJ,INLAYS
!
       NOD=0
       COUNT=0
       INLAYS=0
       DO ILAYER=1,NLAYERS
         DO II=1,NLAY(ILAYER)
           NOD=NOD+1
           FINDRM(NOD)=COUNT+1
           IF(ILAYER.EQ.1) THEN
           ELSE
             DO JJ=1,NLAY(ILAYER-1)
               COUNT=COUNT+1
               COLM(COUNT)=INLAYS-NLAY(ILAYER-1)+JJ
!               write(*,*) 'COUNT,COLM(COUNT):',COUNT,COLM(COUNT)
!               write(*,*) 'INLAYS-NLAY(ILAYER-1)+JJ:', &
!                        INLAYS,NLAY(ILAYER-1),JJ
!               STOP 7
             END DO
           ENDIF
         END DO
         INLAYS=INLAYS+NLAY(ILAYER) 
       END DO
       IF(NONODS.NE.NOD) THEN
          write(*,*) 'NONODS WAS NOT SET TO CORRECT VALUE'
          write(*,*) 'NOD,NONODS:',NOD,NONODS
          STOP 77
       ENDIF
       NCOLM=COUNT
       FINDRM(NONODS+1)=NCOLM+1
!       write(*,*) 'NONODS,NCOLM=',NONODS,NCOLM
!        stop 4
!
!       DO NOD=1,NONODS
!         write(*,*) 'NOD,FINDRM(NOD),FINDRM(NOD+1)-1:', &
!                  NOD,FINDRM(NOD),FINDRM(NOD+1)-1
!         write(*,*) 'COLM:', &
!                (COLM(COUNT),COUNT=FINDRM(NOD),FINDRM(NOD+1)-1)
!       END DO
!       STOP 33
!
       DO COUNT=1,NCOLM
         WEIGHT(COUNT)=0.
       END DO
       RETURN
       END SUBROUTINE SPAWEI
!
!
!
!
       SUBROUTINE GETNEUVALS(NEUVAL,WEIGHT,NONODS, MXNCOLM, &
               NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If NOREXP=1 the dont find exponent of output. 
       IMPLICIT NONE
       LOGICAL OUTEXP
!       PARAMETER(OUTEXP=.FALSE.)
!       PARAMETER(OUTEXP=.true.)
       INTEGER :: NLAYER1,NLAFIN,NONODS,NCOLM, MXNCOLM
       REAL :: NEUVAL(NONODS),WEIGHT(NCOLM)
       INTEGER :: FINDRM(NONODS+1),COLM(MXNCOLM)
       INTEGER :: NOREXP
! LOCAL VARIABLES...
       REAL :: SUM
       INTEGER :: NOD,COUNT
       OUTEXP=(NOREXP.EQ.0)
       DO NOD=NLAYER1+1,NONODS
         SUM=0.
         DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
           SUM=SUM+WEIGHT(COUNT)*NEUVAL(COLM(COUNT))
         END DO
         IF(NOD.EQ.NONODS) THEN
           IF(OUTEXP) THEN
              NEUVAL(NOD)=1./(1.+EXP(-SUM))
           ELSE
              NEUVAL(NOD)=SUM
           ENDIF
         ELSE
           NEUVAL(NOD)=1./(1.+EXP(-SUM))
         ENDIF
       END DO
       RETURN
       END SUBROUTINE GETNEUVALS
!
!
!
!
       SUBROUTINE EOS_NEU(TEMP,PRES,DENSTY, &
           MINT,MAXT,MINP,MAXP,MIND,MAXD, &
           WEIGHT,NONODS, MXNCOLM, &
           NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
! t=(max-min)*scaled+min
! This is the neural network EoS calculates density DENSTY 
! given the temp TEMP and pressure PRES. 
! MINT=min value of temp in data base used for training.
       IMPLICIT NONE
       INTEGER :: MXNODS
       PARAMETER(MXNODS=500)
       REAL :: NEUVAL(MXNODS)  
       INTEGER :: NLAYER1,NLAFIN,NONODS,NCOLM, MXNCOLM
       INTEGER :: FINDRM(NONODS+1),COLM(MXNCOLM)
       INTEGER :: NOREXP
       REAL :: TEMP,PRES,DENSTY
       REAL :: MINT,MAXT,MINP,MAXP,MIND,MAXD
       REAL :: WEIGHT(NCOLM)
! Local variables...
       REAL :: SCALT,SCALP,SCALD
         SCALT=(TEMP-MINT)/(MAXT-MINT)  
         SCALP=(PRES-MINP)/(MAXP-MINP) 
         NEUVAL(1)=SCALT
         NEUVAL(2)=SCALP
         CALL GETNEUVALS(NEUVAL,WEIGHT,NONODS, MXNCOLM, &
               NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
         SCALD=NEUVAL(NONODS)
! Unscale density SCALD to produce DENSTY
         DENSTY=(MAXD-MIND)*SCALD+MIND
         RETURN
         END SUBROUTINE EOS_NEU
!
!
     
        

!!!
!!! In this file, there are all functions and subroutines used to calculate density
!!! at a given temperature and pressure in the following order:
!!! 
!!!     FUNCTION              ARGUMENTS                          RETURNS
!!!========================================================================================
!!!========================================================================================
!!!
!!!1.   TSATUO2               P [g/(cm.s^2)]               Melting temperature for UO2/MOX
!!!                         (deviation in pressure)        after a perturbation in pressure
!!!
!!!----------------------------------------------------------------------------------------
!!!
!!!2.   TSATBOILUO2           P [g/(cm.s^2)]               Boiling temperature for UO2/MOX
!!!                         (deviation in pressure)        after a perturbation in pressure
!!!
!!!----------------------------------------------------------------------------------------
!!!
!!!3.   TSATZNST             P, PATMOS [g/(cm.s^2)],       Saturation temperature for ZnSt
!!!                          TSAT [K]                      as a function of the pressure
!!!                         (deviation in pressure,
!!!                          atmospheric pressure,
!!!                          saturation temperature
!!!                            at normal conditions)
!!!
!!!-----------------------------------------------------------------------------------------
!!!
!!!4.   DENSTYSOLMOX         P [g/(cm.s^2)], T [K],        Density of solid UO2/MOX
!!!                          FRAC
!!!                         (deviation in pressure, 
!!!                         temperature and amount of
!!!                         PuO2 in the solid solution
!!!
!!!-----------------------------------------------------------------------------------------
!!!
!!!5.   DENSTYLIQMOX         P [g/(cm.s^2)], T [K]        Density of liquid UO2/MOX
!!!                         (deviation in pressure and 
!!!                          temperature)
!!!
!!!------------------------------------------------------------------------------------------
!!!
!!!6.   DENSTYVAPMOX         P [g/(cm.s^2)], T [K]        Density of vapor UO2/MOX
!!!                         (deviation in pressure and 
!!!                          temperature)
!!!
!!!------------------------------------------------------------------------------------------
!!!
!!!7.   DENSTYSOLZNST        P [g/(cm.s^2)], T [K]        Density of solid ZnSt
!!!                         (deviation in pressure and 
!!!                          temperature)
!!!
!!!------------------------------------------------------------------------------------------
!!!
!!!8.   DENSTYLIQZNST       P [g/(cm.s^2)], T [K],        Density of liquid ZnSt
!!!                         TMELT [K]
!!!                         (deviation in pressure, 
!!!                          temperature and melting
!!!                          temperature)
!!!
!!!------------------------------------------------------------------------------------------
!!!
!!!9.   DENSTYVAPZNST       P [g/(cm.s^2)], T [K],        Density of vapor ZnSt
!!!                         PATMOS [g/(cm.s^2)]
!!!                         (deviation in pressure, 
!!!                          temperature and atmosphe
!!!                          ric pressure)
!!!
!!!------------------------------------------------------------------------------------------
!!!
!!!10.  SOUNDDEP            P [g/(cm.s^2)] ,T [K],         Speed of sound at current 
!!!                         TMELT, TBOIL, SOUND [cm/s]     temperature and pressure 
!!!                        (deviation in pressure,         conditions
!!!                        current, melting and 
!!!                        boiling temperature and speed 
!!!                       of sound at normal conditions) 
!!!
!!!  ===========================> SUBROUTINE: STIRREDMOX <====================================
!!!
!!! This subrt. calculates the density of MOX (solution of UO2 and PuO2) and ZnSt in a wide
!!!  range of temperature and pressure. The arguments are:
!!! 
!!! 1. P: deviation of atmospheric pressure (PATMOS) in  g/(cm.s^2)
!!! 2. PATMOS: atmospheric pressure (PATMOS) in  g/(cm.s^2)
!!! 3. T: temperature in Kelvin 
!!! 4. EQNOPT: integer that defines if the subrt will calculate density of MOX (=7) or
!!!            ZnSt (=8)
!!! 5. FRAC: mole fraction of PuO2 in the solid MOX solution
!!!        0. <= FRAC <= 1., if there is NO PuO2, then FRAC=0., else if there is
!!!        ONLY PuO2 (NO UO2), then FRAC=1.
!!!
!!!  ======>>>>>>>>>>>  The subroutine will return:  <<<<<<<<<<<<<<<<=======
!!! 
!!! 1. PCOEF
!!! 2. SCOEF
!!!
!!! As so as density will be calculated as
!!!
!!!                 DENSTY = PCOEF * (P + PATMOS) + SCOEF
!!!========================================================================================


REAL FUNCTION TSATUO2(P)
  ! This function calculates the new melting point of UO2 as a function of the 
  ! pressure as reported by 
  ! Musella et al. (1999), 'Dependence of the Melting Temperature on Pressure up 
  ! to 2000 bar in Uranium Dioxide, Tungsten and Graphite', Int. J. Thermophysics,
  ! 20:1177-1188
  ! P: the deviation in pressure of the atmospheric pressure [g/(cm.s^2)]
  ! T: K
  IMPLICIT NONE
  REAL::P
  REAL,PARAMETER::TMELT=3120.,ALPHA=0.04157E-6 ! Slope for the Delta_T=alpha*Delta_P
  !
  TSATUO2=TMELT+P*ALPHA
  !
  RETURN
END FUNCTION TSATUO2
!
!==============================================================================
!

REAL FUNCTION TSATBOILUO2(P)
  ! This function calculates the new boiling point of UO2 as a function of the 
  ! pressure as reported by 
  ! Musella et al. (1999), 'Dependence of the Melting Temperature on Pressure up 
  ! to 2000 bar in Uranium Dioxide, Tungsten and Graphite', Int. J. Thermophysics,
  ! 20:1177-1188
  ! for the melting point
  ! P: the deviation in pressure of the atmospheric pressure [g/(cm.s^2)]
  ! T: K
  IMPLICIT NONE
  REAL::P
  REAL,PARAMETER::TBOIL=3815.,ALPHA=0.04157E-6 ! Slope for the Delta_T=alpha*Delta_P
  !
  TSATBOILUO2=TBOIL+P*ALPHA
  !
  RETURN
END FUNCTION TSATBOILUO2
!
!==============================================================================
!

REAL FUNCTION TSATZNST(P,PATMOS,TSAT)
  ! This function calculates the new melting point of ZnSt as a function of the 
  ! pressure 
  ! P: the deviation in pressure of the atmospheric pressure [g/(cm.s^2)]
  ! PATMOS: atmospheric pressure
  ! T: K
  IMPLICIT NONE
  REAL::P,PATMOS,TSAT
  REAL::AUX
  !
  AUX=7./3.
  TSATZNST=TSAT-P/(PATMOS**AUX)
  !
  RETURN
END FUNCTION TSATZNST
!
!==============================================================================
!

REAL FUNCTION DENSTYSOLMOX(P,T,FRAC)
  ! This function calculates the density of solid UO2 as described by
  ! Carbajo et al. (2001), 'A review of the thermophysical properties of MOX and
  ! UO2 fuels', J. Nuclear Materials, 299:181-198
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  ! FRAC is the amount of PuO2 in percentage, i.e., 
  !    0.0 <= FRAC <= 1.0
  IMPLICIT NONE
  REAL::P,PATMOS,T,FRAC
  !
  REAL,PARAMETER::KAPPA=1.0E-3
  ! Variables for MOX (UO2 and PuO2)
  REAL,PARAMETER::U11=9.9734E-1,U12=9.8020E-6,U13=-2.7050E-10,U14=4.2910E-13
  REAL,PARAMETER::U21=9.9672E-1,U22=1.1790E-5,U23=-2.4290E-9,U24=1.219E-12
  REAL,PARAMETER::DENSTUO2=10.963 ! g/cm3
  !  REAL,PARAMETER::TMOXMELT0=3120.0,TMOXVAP0=3815.
  ! External functions - no, part of module, so don't declare here!
!  REAL::PHIT,SIGNMOX
  !
  REAL::DSTY0,DSTY
  REAL::AUX
  !
  ! Correcting the density of MOX at 273 K
  !
  DSTY0=DENSTUO2+0.497*FRAC
  !
  ! Calculating the density
  !
  IF((T>=273.).AND.(T<923.15))THEN
     DSTY=DSTY0*PHIT(U11,U12,U13,U14,T)
  ELSEIF(T>=923.15)THEN
     DSTY=DSTY0*PHIT(U21,U22,U23,U24,T)
  ENDIF
  AUX=1./11.
  DENSTYSOLMOX=DSTY+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)

  RETURN
END FUNCTION DENSTYSOLMOX
!
!========================================================================

REAL FUNCTION DENSTYLIQMOX(P,T)
  ! This function calculates the density of liquid UO2 as described by
  ! Fink (2000), 'Thermophysical properties of uranium dioxide', J. Nuclear
  ! Materials, 279:1:18
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  IMPLICIT NONE
  REAL::P,T
  REAL,PARAMETER::TMOXMELT0=3120.,KAPPA=1.0E-5
  ! External function - no, part of module, so don't declare here!
!  REAL::SIGNMOX
  !
  REAL::DSTY,AUX
  !
  AUX=1./11.
  DSTY=8860.-0.9285*(T-TSATUO2(P))
  DSTY=DSTY/1000.
  DENSTYLIQMOX=DSTY+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)
  !DENSTYLIQMOX=DENSTYLIQMOX/1000.
  RETURN
END FUNCTION DENSTYLIQMOX
!
!========================================================================
!

REAL FUNCTION DENSTYVAPMOX(P,T)
  ! This function calculates the density of vapor UO2 as described by
  ! Iosilevski et al., Int. J. Thermophysics, 22:1253-1264(2001)
  ! In this case, Iosilevski reported a set of values of density and
  ! we turned it into a second order polynomium
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  IMPLICIT NONE
  REAL::P,T
  REAL,PARAMETER::TMOXMELT=3120.,KAPPA=1.0E-6
  REAL,PARAMETER::A=5.e-6,B=-6.7055225E-07,C=2.0505806E-08
  REAL::DT,AUX!,SIGNMOX
  !
  DT=ABS(T-TMOXMELT)
  DENSTYVAPMOX= A + B*DT + C*(DT**2)
  AUX=1./17.
  DENSTYVAPMOX=DENSTYVAPMOX+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)
  !
  RETURN
END FUNCTION DENSTYVAPMOX
!
!====================================================================
!

REAL FUNCTION DENSTYSOLZNST(P,T)
  ! This function calculates the density of solid ZnSt, in this case,
  ! as we did not find any reliable EOS for plasticizers, such as the 
  ! ZnSt, we decided to use a linear function that would fit some 
  ! density data available in the literature
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  IMPLICIT NONE
  REAL::T,P
  REAL,PARAMETER::DENSTYZNST=1.10,KAPPA=1.0E-3,TCNTP=298.15
  REAL::AUX,DSTY!,SIGNMOX
  !
  AUX=1./6.
  DSTY=DENSTYZNST-1.E-5*(T-TCNTP)
  DENSTYSOLZNST=DSTY+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)
  ! 
  RETURN
END FUNCTION DENSTYSOLZNST
!
!====================================================================
!

REAL FUNCTION DENSTYLIQZNST(P,T,TMELT)
  ! This function calculates the density of liquid ZnSt, in this case,
  ! as we did not find any reliable EOS for plasticizers, such as the 
  ! ZnSt, we decided to use a linear function that would fit some 
  ! density data available in the literature
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  ! TMELT is the melting temperature of ZnSt
  IMPLICIT NONE
  REAL::P,T
  REAL::TMELT
  REAL,PARAMETER::DENSTYZNST=0.6,KAPPA=1.E-3!,TMELT=383.15
  REAL::AUX,DSTY!,SIGNMOX
  !
  AUX=1./6.
  DSTY=DENSTYZNST-1.E-5*(T-TMELT)
  DENSTYLIQZNST=DSTY+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)
  !
  RETURN
END FUNCTION DENSTYLIQZNST
!
!========================================================================
!

REAL FUNCTION DENSTYVAPZNST(P,T,PATMOS)
  ! This function calculates the density of vapor ZnSt, in this case,
  ! as we did not find any reliable EOS for plasticizers, such as the 
  ! ZnSt, we decided to use the ideal gas law. 
  ! The resulting density is corrected for the used pressure with a small
  ! perturbation.
  ! P is the deviation of the atmospheric pressure in g/(cm.s^2)
  ! T is the temperature in Kelvin
  ! PATOMS is the atmospheric pressure
  ! GASCON is gas constant divided by the molecular weight
  IMPLICIT NONE
  REAL::P,T,PATMOS
  REAL,PARAMETER::GASCON=131488.31,KAPPA=1.E-3
  REAL::AUX,DSTY!,SIGNMOX
  !
  AUX=1./20.
  DSTY=PATMOS/(GASCON*T)
  DENSTYVAPZNST=DSTY+KAPPA*SIGNMOX(P)*(ABS(P)**AUX)!+KAPPA*(P**AUX)
  !
  RETURN
END FUNCTION DENSTYVAPZNST
!
!========================================================================
!

REAL FUNCTION SOUNDDEP(P,T,TMELT,TBOIL,SOUND)
  ! This external function corrects the speed of sound with a linear
  ! function for temperature and pressure conditions.  Physical behaviour
  ! was obtained from
  ! 1. Bohn D. A. (1988) 'Environmental Effects on the Speed of Sound",
  !    J. Audio Eng. Soc., Vol 36, No. 4.
  ! 2. Willard G. W. (1947) 'Temperature Coefficient of Ultrasonic Velocity
  !    in Solutions', J. Acoust. Soc. Amer. 19:235-241.
  !
  IMPLICIT NONE
  REAL::P,T,TMELT,TBOIL,SOUND
  REAL,PARAMETER::A=0.0281,B=0.281
  REAL,PARAMETER::T0=298.15
  REAL::SP
  !
  IF(T<TMELT)THEN
     SOUNDDEP=SOUND+B*P
  ELSEIF((T>=TMELT).AND.(T<TBOIL))THEN
     SOUNDDEP=SOUND-A*(T-TMELT)+B*P
  ELSE
     SP=SOUND-A*(TBOIL-TMELT)+B*P
     SOUNDDEP=SP+A*0.01*(T-TBOIL)+B*P/10.
  ENDIF
  !
  RETURN
END FUNCTION SOUNDDEP
!
!========================================================================
!

REAL FUNCTION PHIT(A,B,C,D,T)
  ! This function calculates a third order polynomium
  ! H(T) = A + B*T + C*T^2 + D*T^3 
  IMPLICIT NONE
  REAL::A,B,C,D,T
  !
  PHIT=A+B*T+C*T**2+D*T**3
  PHIT=1./PHIT**3
  RETURN
END FUNCTION PHIT
!
!====================================================================
!

REAL FUNCTION SIGMOIDL(Y,Y0,WIDTH,LOWMAG,UPRMAG)
  !
  ! The sigmoid function, varies between (LOWMAG,UPRMAG).
  ! Y is the variable of the function, Y0 is the centre
  ! of the function and WIDTH is the width of the sigmoid
  ! function.  The function looks like:
  !
  !       --------
  !               /
  !                /
  !                --------
  !
  IMPLICIT NONE
  REAL::Y,Y0,WIDTH,LOWMAG,UPRMAG
  REAL::ALPHA!,EXPREP
  !
  IF(Y-Y0<-0.9*WIDTH)THEN
     SIGMOIDL=UPRMAG
  ELSEIF(Y-Y0>0.9*WIDTH)THEN
     SIGMOIDL=LOWMAG
  ELSE
     ALPHA=35.0/(WIDTH)
     SIGMOIDL=(UPRMAG-LOWMAG)/(1.0+EXP(ALPHA*(Y-Y0)))+LOWMAG
  ENDIF
  !
  RETURN
  !
END FUNCTION SIGMOIDL
!
!===============================================================================
!
!
!=====================================================================================
!

REAL FUNCTION SIGNMOX(A)
  ! This external function returns either -1. or +1. accordingly with the real
  ! variable A, signal
  !
  IMPLICIT NONE
  REAL::A
  !
  SIGNMOX=1.
  IF(A<0.)SIGNMOX=-1.
  !
  RETURN
END FUNCTION SIGNMOX
!
!
!
!

      SUBROUTINE COMBSIMPLEMOXZNST(EQNSTA,PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS)
      IMPLICIT NONE
      INTEGER EQNSTA
      REAL PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS
! this sub calculates the simplified EOS for ZNST and MOX 
! The options are: 
!      IF(OPT.EQ.1) THEN
! Constant density
! options EQNSTA: 14(znst) or 24(mox)
!      ELSE IF(OPT.EQ.2) THEN
! incompressible but density changes with temp (no vapourisation)...
! options 15(znst) or 25(mox)
!
!      ELSE IF(OPT.EQ.3) THEN
! incompressible but density changes with temp (WITH VAPOUR)...
! options 16(znst) or 26(mox)

!      ELSE IF(OPT.EQ.4) THEN
! incompressible but density changes with temp (with ideal gas vapur)...
! options 17(znst) or 27(mox)

!      ELSE IF(OPT.EQ.5) THEN
! incompressible but density changes with temp (with compressible vapur)...
! has a max. gas density of DENLIQ
! options 18(znst) or 28(mox)
! Local variables...
      INTEGER OPT
      IF((EQNSTA.GE.14).AND.(EQNSTA.LE.23)) THEN
! ZNST EOS...
          OPT=EQNSTA-13
          CALL SIMPLEZNST(OPT,PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS)
      ELSE
! MOX EOS...
          OPT=EQNSTA-23
          CALL SIMPLEMOX(OPT,PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS)
      ENDIF
      RETURN
      END SUBROUTINE COMBSIMPLEMOXZNST
!
!
!
      SUBROUTINE SIMPLEMOX(OPT,PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS)
      IMPLICIT NONE
      INTEGER OPT
      REAL DENLIQ,DENSOL,GASCON
      PARAMETER(DENLIQ=8.89,DENSOL=10.9564,GASCON=2478.505)
! Melting and Boiling temperature in Kelvin
! Variables used for MOX calculations
          REAL DELTAMOX,TMELTMOX,TBOILMOX
          PARAMETER(DELTAMOX=10.0)
          PARAMETER(TMELTMOX=3120.0,TBOILMOX=3815.0)
      REAL PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS
! local variables...
      REAL DENGAS,DENGAS2
!  the sigmoid function looks like:
!
!       --------
!               /
!                /
!                --------
      REAL G1,G2
!      G1=SIGMOIDL(TGNEWK,TMELTMOX,5.*DELTAMOX,0.0,1.0)
!      G2=SIGMOIDL(TGNEWK,TBOILMOX,5.*DELTAMOX,0.0,1.0)
      G1=SIGMOIDL(TGNEWK,TMELTMOX,10.*DELTAMOX,0.0,1.0)
      G2=SIGMOIDL(TGNEWK,TBOILMOX,10.*DELTAMOX,0.0,1.0)
      
      CALL SIMPLEMOXZNST(OPT,G1,G2,DENLIQ,DENSOL,GASCON, &
                         PNEW,PATMOS,TGNEWK,PCOEF,DENS)
      
        SCOEF=DENS-PCOEF*(PNEW+PATMOS) 
      RETURN
      END SUBROUTINE SIMPLEMOX
!
!
!
!
      SUBROUTINE SIMPLEZNST(OPT,PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS)
      IMPLICIT NONE
      INTEGER OPT
      REAL DENLIQ,DENSOL,GASCON
!      PARAMETER(DENLIQ=0.6,DENSOL=1.1,GASCON=93354.83)
!      PARAMETER(DENLIQ=0.7,DENSOL=1.1,GASCON=93354.83)
! use 0.8...
      PARAMETER(DENLIQ=0.8,DENSOL=1.1,GASCON=93354.83)
!      PARAMETER(DENLIQ=0.9,DENSOL=1.1,GASCON=93354.83)
!      PARAMETER(DENLIQ=1.1,DENSOL=1.1,GASCON=93354.83)
!      PARAMETER(DENLIQ=1.1,DENSOL=1.1,GASCON=93354.83)
! Melting and Boiling temperature in Kelvin
! Variables used for ZnSt calculations
! Variables used for ZnSt calculations
          REAL TMELTZNST,TBOILZNST,DELTAZNST
          PARAMETER(TMELTZNST=383.15,TBOILZNST=595.15)
          PARAMETER(DELTAZNST=30.0)
      REAL PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS
! local variables...
      REAL DENGAS
!  the sigmoid function looks like:
!
!       --------
!               /
!                /
!                --------
      REAL G1,G2
      G1=SIGMOIDL(TGNEWK,TMELTZNST,2.*DELTAZNST,0.0,1.0)
      G2=SIGMOIDL(TGNEWK,TBOILZNST,2.*DELTAZNST,0.0,1.0)
      
      CALL SIMPLEMOXZNST(OPT,G1,G2,DENLIQ,DENSOL,GASCON, &
                         PNEW,PATMOS,TGNEWK,PCOEF,DENS)
      
        SCOEF=DENS-PCOEF*(PNEW+PATMOS) 
      RETURN
      END SUBROUTINE SIMPLEZNST
!
!
!      
!
      SUBROUTINE SIMPLEMOXZNST(OPT,G1,G2,DENLIQ,DENSOL,GASCON, &
                               PNEW,PATMOS,TGNEWK,PCOEF,DENS)
      IMPLICIT NONE
      INTEGER OPT
      REAL G1,G2
      REAL DENLIQ,DENSOL,GASCON
! Melting and Boiling temperature in Kelvin
! Variables used for MOX calculations
      REAL PNEW,PATMOS,TGNEWK,PCOEF,SCOEF,DENS
!  the sigmoid function looks like (G1,G2):
!
!       --------
!               /
!                /
!                --------
! local variables...
      REAL DENGAS,DENGAS2,Q,Z
      IF(OPT.EQ.1) THEN
! Constant density
! options 14(znst) or 24(mox)
        DENS=DENSOL
        PCOEF=0.0
      ELSE IF(OPT.EQ.2) THEN
! incompressible but density changes with temp (no vapourisation)...
! options 15 or 25
!        DENGAS=(PNEW+PATMOS)/(GASCON*TGNEWK)
        DENGAS=DENLIQ
        DENS=(DENSOL*G1+DENLIQ*(1.-G1))*G2+DENGAS*(1.0-G2)
!        PCOEF=(1.0-G2)/(GASCON*TGNEWK)
        PCOEF=0.0
      ELSE IF(OPT.EQ.3) THEN
! incompressible but density changes with temp (WITH VAPOUR)...
! options 16 or 26
!        DENGAS=(PNEW+PATMOS)/(GASCON*TGNEWK)
        DENGAS=0.1
        DENS=(DENSOL*G1+DENLIQ*(1.-G1))*G2+DENGAS*(1.0-G2)
!        PCOEF=(1.0-G2)/(GASCON*TGNEWK)
        PCOEF=0.0
      ELSE IF(OPT.EQ.4) THEN
! incompressible but density changes with temp (with ideal gas vapur)...
! options 17 or 27
        DENGAS=(PNEW+PATMOS)/(GASCON*TGNEWK)
        DENS=(DENSOL*G1+DENLIQ*(1.-G1))*G2+DENGAS*(1.0-G2)
        PCOEF=(1.0-G2)/(GASCON*TGNEWK)
      ELSE IF(OPT.EQ.5) THEN
! incompressible but density changes with temp (with compressible vapur)...
! has a max. gas density of DENLIQ
! options 18 or 28
        DENGAS2=(PNEW+PATMOS)/(GASCON*TGNEWK)
        DENGAS=(DENLIQ*DENGAS2)/(DENGAS2+DENLIQ)
        DENS=(DENSOL*G1+DENLIQ*(1.-G1))*G2+DENGAS*(1.0-G2)
        Z=DENGAS2
        Q=1./(GASCON*TGNEWK)
        PCOEF=(1.0-G2)*(-DENLIQ*Z/(Z+DENLIQ)**2 &
                        +DENLIQ/(Z+DENLIQ)      )*Q
      ENDIF
      
      SCOEF=DENS-PCOEF*(PNEW+PATMOS) 
      
      RETURN
      END SUBROUTINE SIMPLEMOXZNST
!    
!        
!
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!  MAIN SUBRT TO CALCULATE DENSITY AS A FUNCTION OF PRESSURE AND TEMPERA
!  TURE AND USING A SIGMOID FUNCTION AND THE EXTERNAL FUNCTIONS LISTED 
!  ABOVE
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
SUBROUTINE STIRREDMOX(EQNOPT,FRAC,PRES,PATMOS,TKELV,PCOEF,SCOEF)
  !===================================================================================
  != This subrt. calculates the density of MOX (solution of UO2 and PuO2) and         =
  != ZnSt in a wide range of temperature and pressure                                 =
  != ==> PRES is the deviation of atmospheric pressure (PATMOS) in  g/cm.s2           =
  != ==> EQNOPT: integer that defines if the subrt will calculate density             =
  !=             for ZnSt (EQNOPT=8) or for MOX(=UO2 and PuO2, EQNOPT=7)              =
  != ==> PCOEF and SCOEF are used to calculate pressure as:                           =
  !=       DENST= SCOEF + PCOEF * (PRES + PATMOS)                                     =
  !=     where PCOEF is the sound speed                                               =
  != ==> FRAC is the mole fraction of PuO2 in the solid MOX solution                  =
  !=     0. <= FRAC <= 1., it is adjusted as so as if there is just                   =
  !=     PuO2 (FRAC=1.), then MOX density will be equal to PuO2, else if              =
  !=     there is just UO2, then the resulting density will be the UO2 density.       =
  != ==> TKELV is the temperature in Kelvin                                           =
  !=                                                                                  =
  !=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX=
  != P.S.: Routine tested for tempertures up to 8000 k and pressure 10000 bar         =
  !=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX=
  ! 
  IMPLICIT NONE
  INTEGER::EQNOPT ! Uranium dioxide (EQNOPT=7), zinc stearate (EQNOPT=8)
  REAL::FRAC ! Mole fraction of PuO2 (0. <= FRAC <= 1.)
  REAL::PRES,PATMOS ! Pressure in g.(cm^-1*s^-2)
  REAL::TKELV ! Temperature in Kelvin
  REAL::PCOEF,SCOEF
  REAL::DENST
  !
  REAL::PRESS ! PRESS = PRES + PATMOS
  !
  ! Variables used for MOX (UO2 + PuO2)
  REAL,PARAMETER::DELTAMOX=10.,SOUNDMOX=3155.0E+2 ! Speed of sound in cm/s
  REAL,PARAMETER::TMELTMOX=3120.0,TBOILMOX=3815.  ! Melting and Boiling temperature in Kelvin
  !
  ! Variables used for ZnSt calculations
  REAL,PARAMETER::TMELTZNST=383.15,TBOILZNST=595.15 ! Melting and Boiling temperature in Kelvin
  REAL,PARAMETER::DELTAZNST=30.,SOUNDZNST=3155.E+2  ! Speed of sound in cm/s
  REAL::TMELT,TBOIL,SOUND
  !
  REAL::DENSTY,DENSTUP,DENSTDW
  !
  PRESS=PRES+PATMOS
  !
  IF(EQNOPT==7)THEN ! ======>      EOS for MOX (UO2 + PuO2)
     TMELT=TSATUO2(PRES)
     TBOIL=TSATBOILUO2(PRES)
     !
     IF(TKELV<(TMELT-DELTAMOX))THEN  ! T <= Tmelt - delta
        DENSTY=DENSTYSOLMOX(PRES,TKELV,FRAC)
        ! 
     ELSEIF((TKELV>=(TMELT-DELTAMOX)).AND.(TKELV<(TMELT+DELTAMOX)))THEN ! Tmelt - delta < T <= Tmelt + delta
        DENSTUP=DENSTYSOLMOX(PRES,TMELT-DELTAMOX,FRAC)
        DENSTDW=DENSTYLIQMOX(PRES,TMELT+DELTAMOX)
        DENSTY=SIGMOIDL(TKELV,TMELT,5.*DELTAMOX,DENSTDW,DENSTUP)
        !
     ELSEIF((TKELV>=(TMELT+DELTAMOX)).AND.(TKELV<(TBOIL-DELTAMOX)))THEN ! Tmelt + delta < T < Tvap - delta
        DENSTY=DENSTYLIQMOX(PRES,TKELV)
        !
     ELSEIF((TKELV>=(TBOIL-DELTAMOX)).AND.(TKELV<(TBOIL+DELTAMOX)))THEN ! Tvap - delta <= T <= Tvap + delta
        DENSTUP=DENSTYLIQMOX(PRES,TBOIL-DELTAMOX)
        DENSTDW=DENSTYVAPMOX(PRES,TBOIL+DELTAMOX)
        DENSTY=SIGMOIDL(TKELV,TBOIL,5.*DELTAMOX,DENSTDW,DENSTUP)
        !
     ELSE  ! T > Tvap + delta
!        DENSTY=DENSTYVAPMOX(PRES,TKELV)
        DENSTY=DENSTYVAPMOX(PRES,TBOIL+DELTAMOX)
        !
     ENDIF
     SOUND=SOUNDDEP(PRES,TKELV,TMELT,TBOIL,SOUNDMOX)
     PCOEF=1.0/(SOUND**2.)
     !
  ELSEIF(EQNOPT==8)THEN !  ======>   EOS for ZnSt
     TMELT=TSATZNST(PRES,PATMOS,TMELTZNST) ! Adjusting Tmelt and Tboil for pressure different
     TBOIL=TSATZNST(PRES,PATMOS,TBOILZNST) ! from atmospheric pressure
     !
     IF(TKELV<=(TMELT-DELTAZNST))THEN ! T <= Tmelt - delta
        DENSTY=DENSTYSOLZNST(PRES,TKELV)
        !
     ELSEIF((TKELV>(TMELT-DELTAZNST)).AND.(TKELV<=(TMELT+DELTAZNST)))THEN ! Tmelt - delta < T <= Tmelt + delta
        DENSTUP=DENSTYSOLZNST(PRES,TMELT-DELTAZNST)
        DENSTDW=DENSTYLIQZNST(PRES,TMELT+DELTAZNST,TMELT)
        DENSTY=SIGMOIDL(TKELV,TMELT,2.*DELTAZNST,DENSTDW,DENSTUP)
        !    DENSTY=SIGMOIDL(TKELV,TMELT,2.*DELTAZNST,0.6,DENSTUP)
        !
     ELSEIF((TKELV>(TMELT+DELTAZNST)).AND.(TKELV<=(TBOIL-DELTAZNST)))THEN ! Tmelt + delta < T <= Tboil - delta
        DENSTY=DENSTYLIQZNST(PRES,TKELV,TMELT+DELTAZNST)
        !
     ELSEIF((TKELV>(TBOIL-DELTAZNST)).AND.(TKELV<=(TBOIL+DELTAZNST)))THEN ! Tboil - delta < T <= Tboil + delta
        DENSTUP=DENSTYLIQZNST(PRES,TKELV,TMELT+DELTAZNST)
        DENSTDW=DENSTYVAPZNST(PRES,TBOIL+DELTAZNST,PATMOS)
        DENSTY=SIGMOIDL(TKELV,TBOIL,DELTAZNST,DENSTDW,DENSTUP)
        !
     ELSE ! T > Tboil + delta
        DENSTY=DENSTYVAPZNST(PRES,TKELV,PATMOS)
        !
     ENDIF
     SOUND=SOUNDDEP(PRES,TKELV,TMELT,TBOIL,SOUNDZNST)
     PCOEF=1.0/(SOUND**2.)
  ENDIF
  SCOEF=DENSTY-PCOEF*PRESS


  RETURN
END SUBROUTINE STIRREDMOX


end module steam_nrs_module

!
