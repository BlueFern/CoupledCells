!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   aut.f90 :     Model AUTO-equations file
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----
!
! Evaluates the algebraic equations or ODE right hand side
!
! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters
!
! Values to be returned :
!      F      :   Equation or ODE right hand side values
!
! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

      DOUBLE PRECISION J_PLC
      DOUBLE PRECISION smc_Ca, smc_SR, smc_Vm, smc_w, smc_IP3
      DOUBLE PRECISION Jsmc_Na_Ca, Jsmc_VOCC, Jsmc_Na_K, Jsmc_Cl, Jsmc_K, Jsmc_IP3, Jsmc_SERCA, &
            Jsmc_CICR, Jsmc_Extrusion, Jsmc_Leak, K_activation, Jsmc_IP3_deg
      DOUBLE PRECISION ec_Ca, ec_ER, ec_Vm, ec_IP3
      DOUBLE PRECISION Jec_IP3, Jec_SERCA, Jec_CICR, Jec_Extrusion, Jec_Leak, Jec_IP3_deg, &
            Jec_NSC, Jec_BK_Ca, Jec_SK_Ca, Jec_Ktot, Jec_Residual, Jec_trivial_Ca

      DOUBLE PRECISION Jsmc_hm_Ca, Jsmc_ht_Ca, Jsmc_hm_IP3, Jsmc_ht_IP3, Vmsmc_hm, Vmsmc_ht
      DOUBLE PRECISION Jec_hm_Ca, Jec_ht_Ca, Jec_hm_IP3, Jec_ht_IP3, Vmec_hm, Vmec_ht, tmp

! Constants for homogenically coupled SMCs.
      DOUBLE PRECISION :: Fi = 0.23, Kri = 1.d0, GCai = 0.00129, vCa1 = 100.d0, vCa2 = -24.d0, &
            RCai = 8.50, GNaCai = 0.00316, cNaCai = 0.5, vNaCai = -30.d0, Bi = 2.025, cbi = 1.d0, &
            CICRi = 55.d0, sci = 2.d0, cci = 0.9, Di = 0.24, vdi = -100.d0, Rdi = 250.d0, Li = 0.025, &
            gama = 1970.d0, FNaK = 0.0432, GCli = 0.00134, vCl = -25.d0, GKi = 0.0046, vKi = -94.d0, &
            lambda = 45.d0, cwi = 0.d0, beta = 0.13, vCa3 = -27.0, RKi = 12.d0, ki = 0.1

! Constants for homogenically coupled ECs.
      DOUBLE PRECISION :: Fj = 0.23, Krj = 1.d0, Bj = 0.5, cbj = 1.d0, CICRj = 5.d0, scj = 2.d0, &
            ccj = 0.9, Dj = 0.24, Lj = 0.025, kj = 0.1, Gcatj = 0.00066, ECa = 50.d0, m3cat = -0.18, &
            m4cat = 0.37, J0j = 0.029, Cmj = 25.8, Gtot = 6927.d0, vKj = -80.d0, a1j = 53.3, a2j = 53.3, &
            bj1 = -80.8, c1j = -0.4, m3b = 1.32e-3, m4b = 0.30, m3s = -0.28, m4s = 0.389, GRj = 955.d0, &
            vrestj = -31.10

      J_PLC = PAR(1)

      smc_Ca = U(1)
      smc_SR = U(2)
      smc_Vm = U(3)
      smc_w = U(4)
      smc_IP3 = U(5)

      ec_Ca = U(6)
      ec_ER = U(7)
      ec_Vm = U(8)
      ec_IP3 = U(9)

! SMC fluxes
      Jsmc_Na_Ca = GNaCai * ((smc_Ca / (smc_Ca + cNaCai)) * (smc_Vm - vNaCai));
      Jsmc_VOCC = GCai * (smc_Vm - vCa1) / (1.0 + (exp(-1.0 * (smc_Vm - vCa2) / RCai)));
      Jsmc_Na_K = FNaK;
      Jsmc_Cl = GCli * (smc_Vm - vCl);
      Jsmc_K = GKi * smc_w * (smc_Vm - vKi);
      Jsmc_IP3 = Fi * (smc_IP3**2 / (Kri**2 + smc_IP3**2));
      Jsmc_SERCA = Bi * (smc_Ca**2 / (cbi**2 + smc_Ca**2));
      Jsmc_CICR = CICRi * (smc_SR**2 / (sci**2 + smc_SR**2)) * (smc_Ca**4 / (cci**4 + smc_Ca**4));
      Jsmc_Extrusion = Di * smc_Ca * (1.0 + ((smc_Vm - vdi) / Rdi));
      Jsmc_Leak = Li * smc_SR;
      K_activation = (smc_Ca**2 + cwi)**2 / ((smc_Ca + cwi)**2 + (beta * (exp(-1.0 * (smc_Vm - vCa3) / RKi))));
      Jsmc_IP3_deg = ki * smc_IP3;

! EC fluxes
      Jec_IP3 = Fj * (ec_IP3**2 / (Krj**2 + ec_IP3**2));
      Jec_SERCA = Bj * (ec_Ca**2 / (cbj**2 + ec_Ca**2));
      Jec_CICR = CICRj * (ec_ER**2 / (scj**2 + ec_ER**2)) * (ec_Ca**4 / (ccj**4 + ec_Ca**4));
      Jec_Extrusion = Dj * ec_Ca;
      Jec_Leak = Lj * ec_ER;
      Jec_IP3_deg = kj * ec_IP3;
      Jec_NSC = Gcatj * (ECa - ec_Vm) * 0.5 * (1.0 + tanh((log10(ec_Ca) - m3cat) / m4cat));
      Jec_BK_Ca = (0.4 / 2.0) * (1.0 + tanh(((log10(ec_Ca) - c1j) * (ec_Vm - bj1) - a1j) / &
            (m3b * (ec_Vm + a2j * (log10(ec_Ca) - c1j) - bj1)**2 + m4b)));
      Jec_SK_Ca = (0.6 / 2.0) * (1.0 + tanh((log10(ec_Ca) - m3s) / m4s));
      Jec_Ktot = Gtot * (ec_Vm - vKj) * (Jec_BK_Ca + Jec_SK_Ca);
      Jec_Residual = GRj * (ec_Vm - vrestj);
      Jec_trivial_Ca = J0j;

! SMC Ca
      F(1)= Jsmc_IP3 - Jsmc_SERCA + Jsmc_CICR - Jsmc_Extrusion + Jsmc_Leak- Jsmc_VOCC + Jsmc_Na_Ca + Jsmc_hm_Ca + Jsmc_ht_Ca;
! SMC Ca SR
      F(2)= Jsmc_SERCA - Jsmc_CICR - Jsmc_Leak;
! SMC Vm
      F(3) = gama * (-Jsmc_Na_K - Jsmc_Cl - (2 * Jsmc_VOCC) - Jsmc_Na_Ca - Jsmc_K) + Vmsmc_hm + Vmsmc_ht;
! SMC K
      F(4) = lambda * (K_activation - smc_w);
! SMC IP3
      F(5) = -Jsmc_IP3_deg + Jsmc_hm_IP3 + Jsmc_ht_IP3;
! EC Ca
      F(6) = Jec_IP3 - Jec_SERCA + Jec_CICR - Jec_Extrusion + Jec_Leak + Jec_NSC + Jec_trivial_Ca + Jec_hm_Ca + Jec_ht_Ca;
! EC Ca ER
      F(7) = Jec_SERCA - Jec_CICR - Jec_Leak;
! EC Vm
      F(8) = ((-1 / Cmj) * (Jec_Ktot + Jec_Residual)) + Vmec_hm + Vmec_ht;
! EC IP3
      F(9) = J_PLC - Jec_IP3_deg + Jec_hm_IP3 + Jec_ht_IP3;


      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
!
! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system
!
! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!
! Note : For time- or space-dependent solutions this subroutine has
!        arguments (NDIM,U,PAR,T), where the scalar input parameter T
!        contains the varying time or space variable value.

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

! Initialize the equation parameters
       PAR(1)= 1.0e-3

! Initialize the solution
       U(1)= 1.0e-3
       U(2)= 1.0e-3
       U(3)= 1.0e-3
       U(4)= 1.0e-3
       U(5)= 0.0
       U(6)= 1.0e-3
       U(7)= 1.0e-3
       U(8)= 1.0e-3
       U(9)= 0.0

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----
!
! Boundary Conditions
!
! Input arguments :
!      NDIM   :   Dimension of the ODE system
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NBC    :   Number of boundary conditions
!      U0     :   State variable values at the left boundary
!      U1     :   State variable values at the right boundary
!
! Values to be returned :
!      FB     :   The values of the boundary condition functions
!
! Normally unused Jacobian arguments : IJAC, DBC (see manual)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

!XXX   FB(1)=
!XXX   FB(2)=

      END SUBROUTINE BCND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!     ---------- ----
!
! Integral Conditions
!
! Input arguments :
!      NDIM   :   Dimension of the ODE system
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NINT   :   Number of integral conditions
!      U      :   Value of the vector function U at `time' t
!
! The following input arguments, which are normally not needed,
! correspond to the preceding point on the solution branch
!      UOLD   :   The state vector at 'time' t
!      UDOT   :   Derivative of UOLD with respect to arclength
!      UPOLD  :   Derivative of UOLD with respect to `time'
!
! Normally unused Jacobian arguments : IJAC, DINT
!
! Values to be returned :
!      FI     :   The value of the vector integrand

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

!XXX   FI(1)=

      END SUBROUTINE ICND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
!     ---------- ----
!
! Defines the objective function for algebraic optimization problems
!
! Supplied variables :
!      NDIM   :   Dimension of the state equation
!      U      :   The state vector
!      ICP    :   Indices of the control parameters
!      PAR    :   The vector of control parameters
!
! Values to be returned :
!      FS      :   The value of the objective function
!
! Normally unused Jacobian argument : IJAC, DFDP

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: FS
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

!XXX   FS=

      END SUBROUTINE FOPT
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

!----------------------------------------------------------------------
! NOTE :
! Parameters set in this subroutine should be considered as ``solution
! measures'' and be used for output purposes only.
!
! They should never be used as `true'' continuation parameters.
!
! They may, however, be added as ``over-specified parameters'' in the
! parameter list associated with the AUTO-Constant NICP, in order to
! print their values on the screen and in the ``p.xxx file.
!
! They may also appear in the list associated with AUTO-Constant NUZR.
!
!----------------------------------------------------------------------
! For algebraic problems the argument U is, as usual, the state vector.
! For differential equations the argument U represents the approximate
! solution on the entire interval [0,1]. In this case its values must
! be accessed indirectly by calls to GETP, as illustrated below.
!----------------------------------------------------------------------
!
! Set PAR(2) equal to the L2-norm of U(1)
!XX       PAR(2)=GETP('NRM',1,U)
!
! Set PAR(3) equal to the minimum of U(2)
!XX       PAR(3)=GETP('MIN',2,U)
!
! Set PAR(4) equal to the value of U(2) at the left boundary.
!XX       PAR(4)=GETP('BV0',2,U)
!
! Set PAR(5) equal to the pseudo-arclength step size used.
!XX       PAR(5)=GETP('STP',1,U)
!
!----------------------------------------------------------------------
! The first argument of GETP may be one of the following:
!        'NRM' (L2-norm),     'MAX' (maximum),
!        'INT' (integral),    'BV0 (left boundary value),
!        'MIN' (minimum),     'BV1' (right boundary value).
!
! Also available are
!   'STP' (Pseudo-arclength step size used).
!   'FLD' (`Fold function', which vanishes at folds).
!   'BIF' (`Bifurcation function', which vanishes at singular points).
!   'HBF' (`Hopf function'; which vanishes at Hopf points).
!   'SPB' ( Function which vanishes at secondary periodic bifurcations).
!----------------------------------------------------------------------

      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
