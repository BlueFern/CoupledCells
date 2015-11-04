!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   coupledCells4.f90 :     Coupled Cells Model with 2 SMCs and 2 ECs.
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      MODULE U_INDICES
! Solution vector indices.
      IMPLICIT NONE
      INTEGER, PARAMETER :: U_STRIDE = 9
      INTEGER, PARAMETER :: UId_smc_Ca = 1, UId_smc_SR = 2, UId_smc_Vm = 3, UId_smc_w = 4, UId_smc_IP3 = 5, &
            UId_ec_Ca = 6, UId_ec_ER = 7, UId_ec_Vm = 8, UId_ec_IP3 = 9
      END MODULE

      MODULE FLUX_INDICES
      IMPLICIT NONE
! Indices for SMC fluxes.
     INTEGER, PARAMETER :: SMC_FLX_STRIDE = 12
      INTEGER, PARAMETER :: JId_smc_Na_Ca = 1, JId_smc_VOCC = 2, JId_smc_Na_K = 3, JId_smc_Cl = 4, JId_smc_K = 5, &
            JId_smc_IP3 = 6, JId_smc_SERCA = 7, JId_smc_CICR = 8, JId_smc_Extrusion = 9, JId_smc_Leak = 10, &
            JId_K_activation = 11, JId_smc_IP3_deg = 12
! Indices for EC fluxes.
     INTEGER, PARAMETER :: EC_FLX_STRIDE = 12
      INTEGER, PARAMETER :: JId_ec_IP3 = 1, JId_ec_SERCA = 2, JId_ec_CICR = 3, JId_ec_Extrusion = 4, JId_ec_Leak = 5, &
            JId_ec_IP3_deg = 6, JId_ec_NSC = 7, JId_ec_BK_Ca = 8, JId_ec_SK_Ca = 9, JId_ec_Ktot = 10, &
            JId_ec_Residual = 11, JId_ec_trivial_Ca = 12
      END MODULE

      MODULE PAR_INDICES
      IMPLICIT NONE
! Parameter indices.
      INTEGER, PARAMETER :: PId_J_PLC = 1, PId_Ca_ht = 2, PId_IP3_ht = 3, PId_Vm_ht = 4, PId_Ca_hm_smc = 5, &
            PId_Ca_hm_ec = 6, PId_IP3_hm_smc = 7, PId_IP3_hm_ec = 8, PId_Vm_hm_smc = 9, PId_Vm_hm_ec = 10
      END MODULE

      MODULE COUPL_INDICES
      IMPLICIT NONE
      INTEGER, PARAMETER :: CPL_STRIDE = 12
! Homo coupling indices.
      INTEGER, PARAMETER :: CPL_J_smc_hm_Ca = 1, CPL_J_ec_hm_Ca = 2, CPL_J_smc_hm_IP3 = 3, &
            CPL_J_ec_hm_IP3 = 4, CPL_Vm_smc_hm = 5, CPL_VMm_ec_hm = 6
! Hetero coupling indices.
      INTEGER, PARAMETER:: CPL_J_smc_ht_Ca = 7, CPL_J_ec_ht_Ca = 8, CPL_J_smc_ht_IP3 = 9, &
            CPL_J_ec_ht_IP3 = 10, CPL_Vm_smc_ht = 11, CPL_Vm_ec_ht = 12
      END MODULE

      MODULE CONSTANTS
      IMPLICIT NONE
! Constants for homogenically coupled SMCs.
      DOUBLE PRECISION, PARAMETER :: Fi = 0.23d0, Kri = 1.d0, GCai = 0.00129d0, vCa1 = 100.d0, vCa2 = -24.d0, &
            RCai = 8.5d0, GNaCai = 0.00316d0, cNaCai = 0.5d0, vNaCai = -30.d0, Bi = 2.025d0, cbi = 1.d0, &
            CICRi = 55.d0, sci = 2.d0, cci = 0.9d0, Di = 0.24d0, vdi = -100.d0, Rdi = 250.d0, Li = 0.025d0, &
            gama = 1970.d0, FNaK = 0.0432d0, GCli = 0.00134d0, vCl = -25.d0, GKi = 0.0046d0, vKi = -94.d0, &
            lambda = 45.d0, cwi = 0.d0, beta = 0.13d0, vCa3 = -27.d0, RKi = 12.d0, ki = 0.1d0

! Constants for homogenically coupled ECs.
      DOUBLE PRECISION, PARAMETER :: Fj = 0.23d0, Krj = 1.d0, Bj = 0.5d0, cbj = 1.d0, CICRj = 5.d0, scj = 2.d0, &
            ccj = 0.9d0, Dj = 0.24d0, Lj = 0.025d0, kj = 0.1d0, Gcatj = 0.00066d0, ECa = 50.d0, m3cat = -0.18d0, &
            m4cat = 0.37d0, J0j = 0.029d0, Cmj = 25.8d0, Gtot = 6927.d0, vKj = -80.d0, a1j = 53.3d0, a2j = 53.3d0, &
            bj1 = -80.8d0, c1j = -0.4d0, m3b = 1.32e-3, m4b = 0.3d0, m3s = -0.28d0, m4s = 0.389d0, GRj = 955.d0, &
            vrestj = -31.1d0
      END MODULE

      SUBROUTINE SMC_FLX(U,U_OFFS,SMC)
! Compute SMC fluxes.
      USE U_INDICES
      USE FLUX_INDICES
      USE CONSTANTS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: U(*)
      INTEGER, INTENT(IN) :: U_OFFS
      DOUBLE PRECISION, INTENT(INOUT) :: SMC(*)
      SMC(JId_smc_Na_Ca) = GNaCai * ((U(U_OFFS + UId_smc_Ca) / (U(U_OFFS + UId_smc_Ca) + cNaCai)) * &
            (U(U_OFFS + UId_smc_Vm) - vNaCai));
      SMC(JId_smc_VOCC) = GCai * (U(U_OFFS + UId_smc_Vm) - vCa1) / (1.0 + (exp(-1.0 * (U(U_OFFS + UId_smc_Vm) - vCa2) / RCai)));
      SMC(JId_smc_Na_K) = FNaK;
      SMC(JId_smc_Cl) = GCli * (U(U_OFFS + UId_smc_Vm) - vCl);
      SMC(JId_smc_K) = GKi * U(U_OFFS + UId_smc_w) * (U(U_OFFS + UId_smc_Vm) - vKi);
      SMC(JId_smc_IP3) = Fi * (U(U_OFFS + UId_smc_IP3)**2 / (Kri**2 + U(U_OFFS + UId_smc_IP3)**2));
      SMC(JId_smc_SERCA) = Bi * (U(U_OFFS + UId_smc_Ca)**2 / (cbi**2 + U(U_OFFS + UId_smc_Ca)**2));
      SMC(JId_smc_CICR) = CICRi * (U(U_OFFS + UId_smc_SR)**2 / &
            (sci**2 + U(U_OFFS + UId_smc_SR)**2)) * (U(U_OFFS + UId_smc_Ca)**4 / (cci**4 + U(U_OFFS + UId_smc_Ca)**4));
      SMC(JId_smc_Extrusion) = Di * U(U_OFFS + UId_smc_Ca) * (1.0 + ((U(U_OFFS + UId_smc_Vm) - vdi) / Rdi));
      SMC(JId_smc_Leak) = Li * U(U_OFFS + UId_smc_SR);
      SMC(JId_K_activation) = (U(U_OFFS + UId_smc_Ca)**2 + cwi)**2 / & 
            ((U(U_OFFS + UId_smc_Ca) + cwi)**2 + (beta * (exp(-1.0 * (U(U_OFFS + UId_smc_Vm) - vCa3) / RKi))));
      SMC(JId_smc_IP3_deg) = ki * U(U_OFFS + UId_smc_IP3);
      END SUBROUTINE SMC_FLX

      SUBROUTINE EC_FLX(U,U_OFFS,EC)
! Compute EC fluxes
      USE U_INDICES
      USE FLUX_INDICES
      USE CONSTANTS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: U(*)
      INTEGER, INTENT(IN) :: U_OFFS
      DOUBLE PRECISION, INTENT(INOUT) :: EC(*)
      EC(JId_ec_IP3) = Fj * (U(U_OFFS + UId_ec_IP3)**2 / (Krj**2 + U(U_OFFS + UId_ec_IP3)**2));
      EC(JId_ec_SERCA) = Bj * (U(U_OFFS + UId_ec_Ca)**2 / (cbj**2 + U(U_OFFS + UId_ec_Ca)**2));
      EC(JId_ec_CICR) = CICRj * (U(U_OFFS + UId_ec_ER)**2 / (scj**2 + U(U_OFFS + UId_ec_ER)**2)) * &
            (U(U_OFFS + UId_ec_Ca)**4 / (ccj**4 + U(U_OFFS + UId_ec_Ca)**4));
      EC(JId_ec_Extrusion) = Dj * U(U_OFFS + UId_ec_Ca);
      EC(JId_ec_Leak) = Lj * U(U_OFFS + UId_ec_ER);
      EC(JId_ec_IP3_deg) = kj * U(U_OFFS + UId_ec_IP3);
      EC(JId_ec_NSC) = Gcatj * (ECa - U(U_OFFS + UId_ec_Vm)) * 0.5 * (1.0 + tanh((log10(U(U_OFFS + UId_ec_Ca)) - m3cat) / m4cat));
      EC(JId_ec_BK_Ca) = (0.4 / 2.0) * (1.0 + tanh(((log10(U(U_OFFS + UId_ec_Ca)) - c1j) * (U(U_OFFS + UId_ec_Vm) - bj1) - a1j) / &
            (m3b * (U(U_OFFS + UId_ec_Vm) + a2j * (log10(U(U_OFFS + UId_ec_Ca)) - c1j) - bj1)**2 + m4b)));
      EC(JId_ec_SK_Ca) = (0.6 / 2.0) * (1.0 + tanh((log10(U(U_OFFS + UId_ec_Ca)) - m3s) / m4s));
      EC(JId_ec_Ktot) = Gtot * (U(U_OFFS + UId_ec_Vm) - vKj) * (EC(JId_ec_BK_Ca) + EC(JId_ec_SK_Ca));
      EC(JId_ec_Residual) = GRj * (U(U_OFFS + UId_ec_Vm) - vrestj);
      EC(JId_ec_trivial_Ca) = J0j;
      END SUBROUTINE EC_FLX

      SUBROUTINE SMC_EC_DERIV(SMC,EC,CPL,CPL_OFFS,U,U_OFFS,PAR,F,F_OFFS)
! Compute all state variables.
      USE U_INDICES
      USE PAR_INDICES
      USE FLUX_INDICES
      USE COUPL_INDICES
      USE CONSTANTS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: SMC(*)
      DOUBLE PRECISION, INTENT(IN) :: EC(*)
      DOUBLE PRECISION, INTENT(IN) :: CPL(*)
      INTEGER, INTENT(IN) :: CPL_OFFS
      DOUBLE PRECISION, INTENT(IN) :: U(*)
      INTEGER, INTENT(IN) :: U_OFFS
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(*)
      INTEGER, INTENT(IN) :: F_OFFS
! SMC Ca
      F(F_OFFS + UId_smc_Ca)= SMC(JId_smc_IP3) - SMC(JId_smc_SERCA) + SMC(JId_smc_CICR) - SMC(JId_smc_Extrusion) + &
            SMC(JId_smc_Leak) - SMC(JId_smc_VOCC) + SMC(JId_smc_Na_Ca) + &
            CPL(CPL_OFFS + CPL_J_smc_hm_Ca) + CPL(CPL_OFFS + CPL_J_smc_ht_Ca);
! SMC Ca SR
      F(F_OFFS + UId_smc_SR)= SMC(JId_smc_SERCA) - SMC(JId_smc_CICR) - SMC(JId_smc_Leak);
! SMC Vm
      F(F_OFFS + UId_smc_Vm) = gama * (-SMC(JId_smc_Na_K) - SMC(JId_smc_Cl) - (2.0d0 * SMC(JId_smc_VOCC)) - SMC(JId_smc_Na_Ca) - &
            SMC(JId_smc_K)) + CPL(CPL_OFFS + CPL_Vm_smc_hm) + CPL(CPL_OFFS + CPL_Vm_smc_ht);
! SMC K
      F(F_OFFS + UId_smc_w) = lambda * (SMC(JId_K_activation) - U(U_OFFS + UId_smc_w));
! SMC IP3
      F(F_OFFS + UId_smc_IP3) = - SMC(JId_smc_IP3_deg) + &
            CPL(CPL_OFFS + CPL_J_smc_hm_IP3) + CPL(CPL_OFFS + CPL_J_smc_ht_IP3);
! EC Ca
      F(F_OFFS + UId_ec_Ca) = EC(JId_ec_IP3) - EC(JId_ec_SERCA) + EC(JId_ec_CICR) - EC(JId_ec_Extrusion) + EC(JId_ec_Leak) + &
            EC(JId_ec_NSC) + EC(JId_ec_trivial_Ca) + &
            CPL(CPL_OFFS + CPL_J_ec_hm_Ca) + CPL(CPL_OFFS + CPL_J_ec_ht_Ca);
! EC Ca ER
      F(F_OFFS + UId_ec_ER) = EC(JId_ec_SERCA) - EC(JId_ec_CICR) - EC(JId_ec_Leak);
! EC Vm
      F(F_OFFS + UId_ec_Vm) = ((-1.0d0 / Cmj) * (EC(JId_ec_Ktot) + EC(JId_ec_Residual))) + &
            CPL(CPL_OFFS + CPL_VMm_ec_hm) + CPL(CPL_OFFS + CPL_Vm_ec_ht);
! EC IP3
      F(F_OFFS + UId_ec_IP3) = PAR(PId_J_PLC) - EC(JId_ec_IP3_deg) + &
            CPL(CPL_OFFS + CPL_J_ec_hm_IP3) + CPL(CPL_OFFS + CPL_J_ec_ht_IP3);
      END SUBROUTINE

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

      USE U_INDICES
      USE PAR_INDICES
      USE COUPL_INDICES
      USE FLUX_INDICES
      USE CONSTANTS

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

      INTEGER :: NUM_PAIRS
      DOUBLE PRECISION :: SMC1(SMC_FLX_STRIDE)
      DOUBLE PRECISION :: SMC2(SMC_FLX_STRIDE)
      DOUBLE PRECISION :: EC1(EC_FLX_STRIDE)
      DOUBLE PRECISION :: EC2(EC_FLX_STRIDE)
      DOUBLE PRECISION :: CPL(2*CPL_STRIDE)

! TODO: This can be used for dynamic memory allocation for the SMC, EC & CPL arrays.
! At the moment it's just sitting here and looking idle.
      NUM_PAIRS = NDIM / U_STRIDE

! SMC & EC fluxes.
      CALL SMC_FLX(U,0*U_STRIDE,SMC1)
      CALL SMC_FLX(U,1*U_STRIDE,SMC2)
      CALL EC_FLX(U,0*U_STRIDE,EC1)
      CALL EC_FLX(U,1*U_STRIDE,EC2)

! Homo coupling 1
      CPL(0 + CPL_J_smc_hm_Ca) = -PAR(PId_Ca_hm_smc) * (U(UId_smc_Ca) - U(U_STRIDE + UId_smc_Ca))
      CPL(0 + CPL_J_ec_hm_Ca) = -PAR(PId_Ca_hm_ec) * (U(UId_ec_Ca) - U(U_STRIDE + UId_ec_Ca))
      CPL(0 + CPL_J_smc_hm_IP3) = -PAR(PId_IP3_hm_smc) * (U(UId_smc_IP3) - U(U_STRIDE + UId_smc_IP3))
      CPL(0 + CPL_J_ec_hm_IP3) = -PAR(PId_IP3_hm_ec) * (U(UId_ec_IP3) - U(U_STRIDE + UId_ec_IP3))
      CPL(0 + CPL_Vm_smc_hm) = -PAR(PId_Vm_hm_smc) * (U(UId_smc_Vm) - U(U_STRIDE + UId_smc_Vm))
      CPL(0 + CPL_VMm_ec_hm) = -PAR(PId_Vm_hm_ec) * (U(UId_ec_Vm) - U(U_STRIDE + UId_ec_Vm))

! Homo coupling 2
      CPL(CPL_STRIDE + CPL_J_smc_hm_Ca) = -PAR(PId_Ca_hm_smc) * (U(U_STRIDE + UId_smc_Ca) - U(UId_smc_Ca))
      CPL(CPL_STRIDE + CPL_J_ec_hm_Ca) = -PAR(PId_Ca_hm_ec) * (U(U_STRIDE + UId_ec_Ca) - U(UId_ec_Ca))
      CPL(CPL_STRIDE + CPL_J_smc_hm_IP3) = -PAR(PId_IP3_hm_smc) * (U(U_STRIDE + UId_smc_IP3) - U(UId_smc_IP3))
      CPL(CPL_STRIDE + CPL_J_ec_hm_IP3) = -PAR(PId_IP3_hm_ec) * (U(U_STRIDE + UId_ec_IP3) - U(UId_ec_IP3))
      CPL(CPL_STRIDE + CPL_Vm_smc_hm) = -PAR(PId_Vm_hm_smc) * (U(U_STRIDE + UId_smc_Vm) - U(UId_smc_Vm))
      CPL(CPL_STRIDE + CPL_VMm_ec_hm) = -PAR(PId_Vm_hm_ec) * (U(U_STRIDE + UId_ec_Vm) - U(UId_ec_Vm))

! Hetero coupling 1
      CPL(0 + CPL_J_smc_ht_Ca) = -PAR(PId_Ca_ht) * &
            (2 * U(UId_smc_Ca) - U(UId_ec_Ca) - U(U_STRIDE + UId_ec_Ca))
      CPL(0 + CPL_J_ec_ht_Ca) = -PAR(PId_Ca_ht) * &
            (2 * U(UId_ec_Ca) - U(UId_smc_Ca) - U(U_STRIDE + UId_smc_Ca))
      CPL(0 + CPL_J_smc_ht_IP3) = -PAR(PId_IP3_ht) * &
            (2 * U(UId_smc_IP3) - U(UId_ec_IP3) - U(U_STRIDE + UId_ec_IP3))
      CPL(0 + CPL_J_ec_ht_IP3) = -PAR(PId_IP3_ht) * &
            (2 * U(UId_ec_IP3) - U(UId_smc_IP3) - U(U_STRIDE + UId_smc_IP3))
      CPL(0 + CPL_Vm_smc_ht) = -PAR(PId_Vm_ht) * &
            (2 * U(UId_smc_Vm) - U(UId_ec_Vm) - U(U_STRIDE + UId_ec_Vm))
      CPL(0 + CPL_Vm_ec_ht) = -PAR(PId_Vm_ht) * &
            (2 * U(UId_ec_Vm) - U(UId_smc_Vm) - U(U_STRIDE + UId_smc_Vm))

! Hetero coupling 2
      CPL(CPL_STRIDE + CPL_J_smc_ht_Ca) = -PAR(PId_Ca_ht) * &
            (2 * U(U_STRIDE + UId_smc_Ca) - U(UId_ec_Ca) - U(U_STRIDE + UId_ec_Ca))
      CPL(CPL_STRIDE + CPL_J_ec_ht_Ca) = -PAR(PId_Ca_ht) * &
            (2 * U(U_STRIDE + UId_ec_Ca) - U(UId_smc_Ca) - U(U_STRIDE + UId_smc_Ca))
      CPL(CPL_STRIDE + CPL_J_smc_ht_IP3) = -PAR(PId_IP3_ht) * &
            (2 * U(U_STRIDE + UId_smc_IP3) - U(UId_ec_IP3) - U(U_STRIDE + UId_ec_IP3))
      CPL(CPL_STRIDE + CPL_J_ec_ht_IP3) = -PAR(PId_IP3_ht) * &
            (2 * U(U_STRIDE + UId_ec_IP3) - U(UId_smc_IP3) - U(U_STRIDE + UId_smc_IP3))
      CPL(CPL_STRIDE + CPL_Vm_smc_ht) = -PAR(PId_Vm_ht) * &
            (2 * U(U_STRIDE + UId_smc_Vm) - U(UId_ec_Vm) - U(U_STRIDE + UId_ec_Vm))
      CPL(CPL_STRIDE + CPL_Vm_ec_ht) = -PAR(PId_Vm_ht) * &
            (2 * U(U_STRIDE + UId_ec_Vm) - U(UId_smc_Vm) - U(U_STRIDE + UId_smc_Vm))

! SMC & EC derivatives.
      CALL SMC_EC_DERIV(SMC1,EC1,CPL,0*CPL_STRIDE,U,0*U_STRIDE,PAR,F,0*U_STRIDE)
      CALL SMC_EC_DERIV(SMC2,EC2,CPL,1*CPL_STRIDE,U,1*U_STRIDE,PAR,F,1*U_STRIDE)

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
      USE PAR_INDICES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

! JPLC start.
       PAR(PId_J_PLC) = 0.d0

! Hetero coupling CASE 1.
       PAR(PId_Ca_ht) = 0.d0
       PAR(PId_IP3_ht) = 0.05d0
       PAR(PId_Vm_ht) = 50.d0

! Homo coupling ZERO.
!       PAR(PId_Ca_hm_smc) = 0.0d0
!       PAR(PId_Ca_hm_ec) = 0.0d0
!       PAR(PId_IP3_hm_smc) = 0.0d0
!       PAR(PId_IP3_hm_ec) = 0.0d0
!       PAR(PId_Vm_hm_smc) = 0.0d0
!       PAR(PId_Vm_hm_ec) = 0.0d0

! Homo coupling CASE 1.
       PAR(PId_Ca_hm_smc) = 0.05d0
       PAR(PId_Ca_hm_ec) = 0.05d0
       PAR(PId_IP3_hm_smc) = 0.05d0
       PAR(PId_IP3_hm_ec) = 0.0d0
       PAR(PId_Vm_hm_smc) = 1000.0d0
       PAR(PId_Vm_hm_ec) = 1000.0d0

! Initialize the solution
       U = 1.0e-3

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
      END SUBROUTINE BCND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      END SUBROUTINE ICND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
      END SUBROUTINE FOPT
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
