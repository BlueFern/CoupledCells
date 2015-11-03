!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   coupledCells4.f90 :     Coupled Cells Model with 2 SMCs and 2 ECs
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

! SMC and EC no. 1 fluxes
      DOUBLE PRECISION Jsmc_Na_Ca_1, Jsmc_VOCC_1, Jsmc_Na_K_1, Jsmc_Cl_1, Jsmc_K_1, Jsmc_IP3_1, Jsmc_SERCA_1, &
            Jsmc_CICR_1, Jsmc_Extrusion_1, Jsmc_Leak_1, K_activation_1, Jsmc_IP3_1_deg_1

      DOUBLE PRECISION Jec_IP3_1, Jec_SERCA_1, Jec_CICR_1, Jec_Extrusion_1, Jec_Leak_1, Jec_IP3_1_deg_1, &
            Jec_NSC_1, Jec_BK_Ca_1, Jec_SK_Ca_1, Jec_Ktot_1, Jec_Residual_1, Jec_trivial_Ca_1

      DOUBLE PRECISION smc_Ca_1, smc_SR_1, smc_Vm_1, smc_w_1, smc_IP3_1
      DOUBLE PRECISION ec_Ca_1, ec_ER_1, ec_Vm_1, ec_IP3_1

      DOUBLE PRECISION Jsmc_hm_Ca_1, Jsmc_hm_IP3_1, Vmsmc_hm_1
      DOUBLE PRECISION Jec_hm_Ca_1, Jec_hm_IP3_1, Vmec_hm_1
 
      DOUBLE PRECISION Jsmc_ht_Ca_1, Jec_ht_Ca_1, Jsmc_ht_IP3_1, Jec_ht_IP3_1, Vmsmc_ht_1, Vmec_ht_1

! SMC and EC no. 1 fluxes
      DOUBLE PRECISION Jsmc_Na_Ca_2, Jsmc_VOCC_2, Jsmc_Na_K_2, Jsmc_Cl_2, Jsmc_K_2, Jsmc_IP3_2, Jsmc_SERCA_2, &
            Jsmc_CICR_2, Jsmc_Extrusion_2, Jsmc_Leak_2, K_activation_2, Jsmc_IP3_2_deg_2

      DOUBLE PRECISION Jec_IP3_2, Jec_SERCA_2, Jec_CICR_2, Jec_Extrusion_2, Jec_Leak_2, Jec_IP3_2_deg_2, &
            Jec_NSC_2, Jec_BK_Ca_2, Jec_SK_Ca_2, Jec_Ktot_2, Jec_Residual_2, Jec_trivial_Ca_2

      DOUBLE PRECISION smc_Ca_2, smc_SR_2, smc_Vm_2, smc_w_2, smc_IP3_2
      DOUBLE PRECISION ec_Ca_2, ec_ER_2, ec_Vm_2, ec_IP3_2

      DOUBLE PRECISION Jsmc_hm_Ca_2, Jsmc_hm_IP3_2, Vmsmc_hm_2
      DOUBLE PRECISION Jec_hm_Ca_2, Jec_hm_IP3_2, Vmec_hm_2

      DOUBLE PRECISION Jsmc_ht_Ca_2, Jec_ht_Ca_2, Jsmc_ht_IP3_2, Jec_ht_IP3_2, Vmsmc_ht_2, Vmec_ht_2

! Parameters
      DOUBLE PRECISION J_PLC
      DOUBLE PRECISION Ca_ht, IP3_ht, Vm_ht
      DOUBLE PRECISION Ca_hm_smc, IP3_hm_smc, Vm_hm_smc, Ca_hm_ec, IP3_hm_ec, Vm_hm_ec

      J_PLC = PAR(1)
      Ca_ht = PAR(2)
      IP3_ht = PAR(3)
      Vm_ht = PAR(4)
      Ca_hm_smc = PAR(5)
      IP3_hm_smc = PAR(6)
      Vm_hm_smc = PAR(7)
      Ca_hm_ec = PAR(8)
      IP3_hm_ec = PAR(9)
      Vm_hm_ec = PAR(10)

      smc_Ca_1 = U(1)
      smc_Ca_2 = U(2)
      smc_SR_1 = U(3)
      smc_SR_2 = U(4)
      smc_Vm_1 = U(5)
      smc_Vm_2 = U(6)
      smc_w_1 = U(7)
      smc_w_2 = U(8)
      smc_IP3_1 = U(9)
      smc_IP3_2 = U(10)

      ec_Ca_1 = U(11)
      ec_Ca_2 = U(12)
      ec_ER_1 = U(13)
      ec_ER_2 = U(14)
      ec_Vm_1 = U(15)
      ec_Vm_2 = U(16)
      ec_IP3_1 = U(17)
      ec_IP3_2 = U(18)

! SMC fluxes 1
      Jsmc_Na_Ca_1 = GNaCai * ((smc_Ca_1 / (smc_Ca_1 + cNaCai)) * (smc_Vm_1 - vNaCai));
      Jsmc_VOCC_1 = GCai * (smc_Vm_1 - vCa1) / (1.0 + (exp(-1.0 * (smc_Vm_1 - vCa2) / RCai)));
      Jsmc_Na_K_1 = FNaK;
      Jsmc_Cl_1 = GCli * (smc_Vm_1 - vCl);
      Jsmc_K_1 = GKi * smc_w_1 * (smc_Vm_1 - vKi);
      Jsmc_IP3_1 = Fi * (smc_IP3_1**2 / (Kri**2 + smc_IP3_1**2));
      Jsmc_SERCA_1 = Bi * (smc_Ca_1**2 / (cbi**2 + smc_Ca_1**2));
      Jsmc_CICR_1 = CICRi * (smc_SR_1**2 / (sci**2 + smc_SR_1**2)) * (smc_Ca_1**4 / (cci**4 + smc_Ca_1**4));
      Jsmc_Extrusion_1 = Di * smc_Ca_1 * (1.0 + ((smc_Vm_1 - vdi) / Rdi));
      Jsmc_Leak_1 = Li * smc_SR_1;
      K_activation_1 = (smc_Ca_1**2 + cwi)**2 / ((smc_Ca_1 + cwi)**2 + (beta * (exp(-1.0 * (smc_Vm_1 - vCa3) / RKi))));
      Jsmc_IP3_1_deg_1 = ki * smc_IP3_1;

! EC fluxes 1
      Jec_IP3_1 = Fj * (ec_IP3_1**2 / (Krj**2 + ec_IP3_1**2));
      Jec_SERCA_1 = Bj * (ec_Ca_1**2 / (cbj**2 + ec_Ca_1**2));
      Jec_CICR_1 = CICRj * (ec_ER_1**2 / (scj**2 + ec_ER_1**2)) * (ec_Ca_1**4 / (ccj**4 + ec_Ca_1**4));
      Jec_Extrusion_1 = Dj * ec_Ca_1;
      Jec_Leak_1 = Lj * ec_ER_1;
      Jec_IP3_1_deg_1 = kj * ec_IP3_1;
      Jec_NSC_1 = Gcatj * (ECa - ec_Vm_1) * 0.5 * (1.0 + tanh((log10(ec_Ca_1) - m3cat) / m4cat));
      Jec_BK_Ca_1 = (0.4 / 2.0) * (1.0 + tanh(((log10(ec_Ca_1) - c1j) * (ec_Vm_1 - bj1) - a1j) / &
            (m3b * (ec_Vm_1 + a2j * (log10(ec_Ca_1) - c1j) - bj1)**2 + m4b)));
      Jec_SK_Ca_1 = (0.6 / 2.0) * (1.0 + tanh((log10(ec_Ca_1) - m3s) / m4s));
      Jec_Ktot_1 = Gtot * (ec_Vm_1 - vKj) * (Jec_BK_Ca_1 + Jec_SK_Ca_1);
      Jec_Residual_1 = GRj * (ec_Vm_1 - vrestj);
      Jec_trivial_Ca_1 = J0j;

! Homo coupling 1
      Jsmc_hm_Ca_1 = -Ca_hm_smc * (smc_Ca_1 - smc_Ca_2)
      Jec_hm_Ca_1 = -Ca_hm_ec * (ec_Ca_1 - ec_Ca_2)
      Jsmc_hm_IP3_1 = -IP3_hm_smc * (smc_Ca_1 - smc_Ca_2)
      Jec_hm_IP3_1 = -IP3_hm_ec * (ec_Ca_1 - ec_Ca_2)
      Vmsmc_hm_1 = -Vm_hm_smc * (smc_Ca_1 - smc_Ca_2)
      Vmec_hm_1 = -Vm_hm_ec * (ec_Ca_1 - ec_Ca_2)

! Hetero coupling 1
      Jsmc_ht_Ca_1 = -Ca_ht * (smc_Ca_1 - ec_Ca_1)
      Jec_ht_Ca_1 = -Ca_ht * (ec_Ca_1 - smc_Ca_1)
      Jsmc_ht_IP3_1 = -IP3_ht * (smc_IP3_1 - ec_IP3_1)
      Jec_ht_IP3_1 = -IP3_ht * (ec_IP3_1 - smc_IP3_1)
      Vmsmc_ht_1 = -Vm_ht * (smc_Vm_1 - ec_Vm_1)
      Vmec_ht_1 = -Vm_ht * (ec_Vm_1 - smc_Vm_1)



! SMC fluxes 2
      Jsmc_Na_Ca_2 = GNaCai * ((smc_Ca_2 / (smc_Ca_2 + cNaCai)) * (smc_Vm_2 - vNaCai));
      Jsmc_VOCC_2 = GCai * (smc_Vm_2 - vCa1) / (1.0 + (exp(-1.0 * (smc_Vm_2 - vCa2) / RCai)));
      Jsmc_Na_K_2 = FNaK;
      Jsmc_Cl_2 = GCli * (smc_Vm_2 - vCl);
      Jsmc_K_2 = GKi * smc_w_2 * (smc_Vm_2 - vKi);
      Jsmc_IP3_2 = Fi * (smc_IP3_2**2 / (Kri**2 + smc_IP3_2**2));
      Jsmc_SERCA_2 = Bi * (smc_Ca_2**2 / (cbi**2 + smc_Ca_2**2));
      Jsmc_CICR_2 = CICRi * (smc_SR_2**2 / (sci**2 + smc_SR_2**2)) * (smc_Ca_2**4 / (cci**4 + smc_Ca_2**4));
      Jsmc_Extrusion_2 = Di * smc_Ca_2 * (1.0 + ((smc_Vm_2 - vdi) / Rdi));
      Jsmc_Leak_2 = Li * smc_SR_2;
      K_activation_2 = (smc_Ca_2**2 + cwi)**2 / ((smc_Ca_2 + cwi)**2 + (beta * (exp(-1.0 * (smc_Vm_2 - vCa3) / RKi))));
      Jsmc_IP3_2_deg_2 = ki * smc_IP3_2;

! EC fluxes 2
      Jec_IP3_2 = Fj * (ec_IP3_2**2 / (Krj**2 + ec_IP3_2**2));
      Jec_SERCA_2 = Bj * (ec_Ca_2**2 / (cbj**2 + ec_Ca_2**2));
      Jec_CICR_2 = CICRj * (ec_ER_2**2 / (scj**2 + ec_ER_2**2)) * (ec_Ca_2**4 / (ccj**4 + ec_Ca_2**4));
      Jec_Extrusion_2 = Dj * ec_Ca_2;
      Jec_Leak_2 = Lj * ec_ER_2;
      Jec_IP3_2_deg_2 = kj * ec_IP3_2;
      Jec_NSC_2 = Gcatj * (ECa - ec_Vm_2) * 0.5 * (1.0 + tanh((log10(ec_Ca_2) - m3cat) / m4cat));
      Jec_BK_Ca_2 = (0.4 / 2.0) * (1.0 + tanh(((log10(ec_Ca_2) - c1j) * (ec_Vm_2 - bj1) - a1j) / &
            (m3b * (ec_Vm_2 + a2j * (log10(ec_Ca_2) - c1j) - bj1)**2 + m4b)));
      Jec_SK_Ca_2 = (0.6 / 2.0) * (1.0 + tanh((log10(ec_Ca_2) - m3s) / m4s));
      Jec_Ktot_2 = Gtot * (ec_Vm_2 - vKj) * (Jec_BK_Ca_2 + Jec_SK_Ca_2);
      Jec_Residual_2 = GRj * (ec_Vm_2 - vrestj);
      Jec_trivial_Ca_2 = J0j;

! Homo coupling 2
      Jsmc_hm_Ca_2 = -Ca_hm_smc * (smc_Ca_2 - smc_Ca_1)
      Jec_hm_Ca_2 = -Ca_hm_ec * (ec_Ca_2 - ec_Ca_1)
      Jsmc_hm_IP3_2 = -IP3_hm_smc * (smc_Ca_2 - smc_Ca_1)
      Jec_hm_IP3_2 = -IP3_hm_ec * (ec_Ca_2 - ec_Ca_1)
      Vmsmc_hm_2 = -Vm_hm_smc * (smc_Ca_2 - smc_Ca_1)
      Vmec_hm_2 = -Vm_hm_ec * (ec_Ca_2 - ec_Ca_1)

! Hetero coupling 2
      Jsmc_ht_Ca_2 = -Ca_ht * (smc_Ca_2 - ec_Ca_2)
      Jec_ht_Ca_2 = -Ca_ht * (ec_Ca_2 - smc_Ca_2)
      Jsmc_ht_IP3_2 = -IP3_ht * (smc_IP3_2 - ec_IP3_2)
      Jec_ht_IP3_2 = -IP3_ht * (ec_IP3_2 - smc_IP3_2)
      Vmsmc_ht_2 = -Vm_ht * (smc_Vm_2 - ec_Vm_2)
      Vmec_ht_2 = -Vm_ht * (ec_Vm_2 - smc_Vm_2)




! SMC Ca 1
      F(1)= Jsmc_IP3_1 - Jsmc_SERCA_1 + Jsmc_CICR_1 - Jsmc_Extrusion_1 + Jsmc_Leak_1 - Jsmc_VOCC_1 + Jsmc_Na_Ca_1 + Jsmc_hm_Ca_1 &
           + Jsmc_ht_Ca_1;
! SMC Ca SR 1
      F(3)= Jsmc_SERCA_1 - Jsmc_CICR_1 - Jsmc_Leak_1;
! SMC Vm 1
      F(5) = gama * (-Jsmc_Na_K_1 - Jsmc_Cl_1 - (2 * Jsmc_VOCC_1) - Jsmc_Na_Ca_1 - Jsmc_K_1) + Vmsmc_hm_1 + Vmsmc_ht_1;
! SMC K 1
      F(7) = lambda * (K_activation_1 - smc_w_1);
! SMC IP3 1
      F(9) = -Jsmc_IP3_1_deg_1 + Jsmc_hm_IP3_1 + Jsmc_ht_IP3_1;
! EC Ca 1
      F(11) = Jec_IP3_1 - Jec_SERCA_1 + Jec_CICR_1 - Jec_Extrusion_1 + Jec_Leak_1 + Jec_NSC_1 + Jec_trivial_Ca_1 + Jec_hm_Ca_1 &
            + Jec_ht_Ca_1;
! EC Ca ER 1
      F(13) = Jec_SERCA_1 - Jec_CICR_1 - Jec_Leak_1;
! EC Vm 1
      F(15) = ((-1 / Cmj) * (Jec_Ktot_1 + Jec_Residual_1)) + Vmec_hm_1 + Vmec_ht_1;
! EC IP3 1
      F(17) = J_PLC - Jec_IP3_1_deg_1 + Jec_hm_IP3_1 + Jec_ht_IP3_1;


! SMC Ca 2
      F(2)= Jsmc_IP3_2 - Jsmc_SERCA_2 + Jsmc_CICR_2 - Jsmc_Extrusion_2 + Jsmc_Leak_2 - Jsmc_VOCC_2 + Jsmc_Na_Ca_2 + Jsmc_hm_Ca_2 &
	+ Jsmc_ht_Ca_2;
! SMC Ca SR 2
      F(4)= Jsmc_SERCA_2 - Jsmc_CICR_2 - Jsmc_Leak_2;
! SMC Vm 2
      F(6) = gama * (-Jsmc_Na_K_2 - Jsmc_Cl_2 - (2 * Jsmc_VOCC_2) - Jsmc_Na_Ca_2 - Jsmc_K_2) + Vmsmc_hm_2 + Vmsmc_ht_2;
! SMC K 2
      F(8) = lambda * (K_activation_2 - smc_w_2);
! SMC IP3 2
      F(10) = -Jsmc_IP3_2_deg_2 + Jsmc_hm_IP3_2 + Jsmc_ht_IP3_2;
! EC Ca 2
      F(12) = Jec_IP3_2 - Jec_SERCA_2 + Jec_CICR_2 - Jec_Extrusion_2 + Jec_Leak_2 + Jec_NSC_2 + Jec_trivial_Ca_2 + Jec_hm_Ca_2 &
	+ Jec_ht_Ca_2;
! EC Ca ER 2
      F(14) = Jec_SERCA_2 - Jec_CICR_2 - Jec_Leak_2;
! EC Vm 2
      F(16) = ((-1 / Cmj) * (Jec_Ktot_2 + Jec_Residual_2)) + Vmec_hm_2 + Vmec_ht_2;
! EC IP3 2
      F(18) = J_PLC - Jec_IP3_2_deg_2 + Jec_hm_IP3_2 + Jec_ht_IP3_2;


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
      !J_PLC
      PAR(1) = 0.0d0
      !Ca_ht
      PAR(2) = 0.d0
      !IP3_ht
      PAR(3) = 0.05d0
      !Vm_ht
      PAR(4) = 50.d0
      !Ca_hm_smc
      PAR(5) = 0.05d0
      !IP3_hm_smc
      PAR(6) = 0.05d0
      !Vm_hm_smc
      PAR(7) = 1000.d0
      !Ca_hm_ec
      PAR(8) = 0.05d0
      !IP3_hm_ec
      PAR(9) = 0.d0
      !Vm_hm_ec
      PAR(10) = 1000.d0


! Initialize the solution
       U(1)= 1.0e-3
       U(2)= 1.0e-3
       U(3)= 1.0e-3
       U(4)= 1.0e-3
       U(5)= 1.0e-3
       U(6)= 1.0e-3
       U(7)= 1.0e-3
       U(8)= 1.0e-3
       U(9)= 1.0e-3
       U(10)= 1.0e-3
       U(11)= 1.0e-3
       U(12)= 1.0e-3
       U(13)= 1.0e-3
       U(14)= 1.0e-3
       U(15)= 1.0e-3
       U(16)= 1.0e-3
       U(17)= 1.0e-3
       U(18)= 1.0e-3

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
