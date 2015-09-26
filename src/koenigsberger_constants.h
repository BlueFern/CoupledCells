#ifndef KOENIGSBERGER_CONSTANTS_H_
#define KOENIGSBERGER_CONSTANTS_H_

/*
 * koenigsberger_constants.h
 *
 * FIXED MODEL PARAMETERS FOR THE KOENIGSBERGER MODEL
 *
 */

double
	/* Constants for homogenically coupled SMCs. */
      Fi = 0.23,  Kri = 1.00,  GCai = 0.00129,  vCa1 = 100.00,
	  vCa2 = -24.00,  RCai = 8.50,  GNaCai = 0.00316,  cNaCai = 0.5,
	  vNaCai = -30.00,  Bi = 2.025,  cbi = 1.0,  CICRi = 55.00,
	  sci = 2.0,  cci = 0.9,  Di = 0.24,  vdi = -100.00,
	  Rdi = 250.00,  Li = 0.025, gama = 1970.00,
	  FNaK = 0.0432,  GCli = 0.00134,  vCl = -25.0,
	  GKi = 0.0046,  vKi = -94.00,  lambda = 45.00,
	  cwi = 0.0,  beta = 0.13,  vCa3 = -27.0,
	  RKi = 12.0,  ki = 0.1,
	/* Constants for homogenically coupled ECs. */
	  Fj = 0.23,  Krj = 1.00,  Bj = 0.5,
	  cbj = 1.0,  CICRj = 5.0,  scj = 2.0,
	  ccj = 0.9,  Dj = 0.24,  Lj = 0.025,
	  kj = 0.1,  Gcatj = 0.66 * 1e-3,  ECa = 50.00,
	  m3cat = -0.18,
	  m4cat = 0.37,  J0j = 0.029,  Cmj = 25.8,
	  Gtot = 6927,   vKj = -80.0,  a1j = 53.3,
	  a2j = 53.3,     bj = -80.8,  c1j = -0.4,
	  m3b = 1.32e-3, m4b = 0.30,   m3s = -0.28,
	  m4s = 0.389,   GRj = 955, vrestj = -31.10,
	/*Intracellular calcium buffering*/
	  k6 = 100.00,  k7 = 300.00,  BT = 120.00;

#endif /* KOENIGSBERGER_CONSTANTS_H_ */
