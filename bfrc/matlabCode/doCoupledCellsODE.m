% Run the coupled cells ODE model and plot smc_Ca against time

clear;

JPLC = 0.275;
Vm_ht = 0;
Ca_ht = 0;
IP3_ht = 0.05;

f = @(t,x) CoupledCellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht);
[t x] = ode15s(f, 0:0.01:1000, [0.1,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12]);

smc_Ca = x(:,1); 

figure(100);
plot(t,smc_Ca);
xlabel('t');