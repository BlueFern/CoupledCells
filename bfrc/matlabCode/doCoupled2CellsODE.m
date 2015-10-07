% Run the 2 coupled cells ODE model and plot smc_Ca_1 and smc_Ca_2 against time

clear;
odeopts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.05, 'Vectorized', 1);

JPLC = 0.275;

% Coupling coefficients
Vm_ht = 0;
Ca_ht = 0;
IP3_ht = 0.05;
Vm_hm_smc = 1000;
Vm_hm_ec = 0;
Ca_hm_smc = 0.05;
Ca_hm_ec = 0;
IP3_hm_smc = 0.05;
IP3_hm_ec = 0.05;

f = @(t,x) Coupled2CellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht,Vm_hm_smc,Ca_hm_smc,IP3_hm_smc,Vm_hm_ec,Ca_hm_ec,IP3_hm_ec);
[t x] = ode15s(f, 0:0.01:1000, [0.3,0.2,-60,0.1,0.4,0.12,0.12,-62,0.12,0.4,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

smc_Ca_1 = x(:,1); 
smc_Ca_2 = x(:,10); 

figure(126);
plot(t,smc_Ca_1, t, smc_Ca_2);
xlabel('t');
legend('smc_{Ca1}','smc_{Ca2}');