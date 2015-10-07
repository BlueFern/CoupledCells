function dxdt = CoupledCellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht)

dxdt = zeros(size(x));
smc_Ca = x(1); 
smc_SR = x(2);
smc_Vm = x(3);
smc_w = x(4);
smc_IP3 = x(5);
ec_Ca = x(6);
ec_SR = x(7);
ec_Vm = x(8);
ec_IP3 = x(9);



dxdt(1) = (0.23 * smc_IP3^2) / (1.00^2 + smc_IP3^2) - (2.025 * smc_Ca^2) / (smc_Ca^2 + 1.00^2)  + (55 * (smc_SR^2 * smc_Ca^4)) / ((2^2 + smc_SR^2) * (0.9^4 + smc_Ca^4)) - 0.24 * smc_Ca * (1 + ((smc_Vm - -100) / 250)) + 0.025 * smc_SR - 0.00129 * (smc_Vm - 100) / (1 + ((exp((-1 * (smc_Vm - -24)) / 8.5)))) + 0.00316 * smc_Ca * (smc_Vm - -30) / (smc_Ca + 0.5)  - Ca_ht*(smc_Ca - ec_Ca);
dxdt(2) = (2.025 * smc_Ca^2) / (smc_Ca^2 + 1.00^2) - (55 * (smc_SR^2 * smc_Ca^4)) / ((2^2 + smc_SR^2) * (0.9^4 + smc_Ca^4)) - 0.025 * smc_SR;                                                                                                                                                                                                                                                   
dxdt(3) = 1970 * (-0.0432 - 0.00134 * (smc_Vm - -25) -  (2 * 0.00129 * (smc_Vm - 100) / (1 + ((exp((-1 * (smc_Vm - -24)) / 8.5))))) - 0.00316 * smc_Ca * (smc_Vm - -30) / (smc_Ca + 0.5) - 0.0046 * smc_w * (smc_Vm - -94))  - Vm_ht*(smc_Vm - ec_Vm);                                                                                                                                           
dxdt(4) = 45 * ((smc_Ca + 0)^2 / ((smc_Ca + 0)^2 + (0.13 * (exp(-1 * (smc_Vm - -27) / 12)))) - smc_w);                                                                                                                                                                                                                                                                                            
dxdt(5) = -0.1 * smc_IP3 - IP3_ht*(smc_IP3 - ec_IP3);                                                                                                                                                                                                                                                                                                                                           
dxdt(6) = (0.23 * ec_IP3^2) / (1^2 + ec_IP3^2) - (0.5 * ec_Ca^2) / (ec_Ca^2 + 1^2) + (5 * ec_SR^2 * ec_Ca^4) / ((ec_SR^2 + 2^2)*(ec_Ca^4 + 0.9^4))  - 0.24 * ec_Ca  + 0.025 * ec_SR  + (0.66*1e-3 * (50 - ec_Vm) * 0.5) * (1 + ((tanh((((log10(ec_Ca)) - -0.18) / 0.37)))))+ 0.029 - Ca_ht*(ec_Ca - smc_Ca);                                                                                      
dxdt(7) = (0.5 * ec_Ca^2) / (ec_Ca^2 + 1^2) - (5 * ec_SR^2 * ec_Ca^4) / ((ec_SR^2 + 2^2)*(ec_Ca^4+0.9^4)) - 0.025 * ec_SR;                                                                                                                                                                                                                                                                        
dxdt(8) = ((-1 / 25.8) * (6927 * (ec_Vm - -80) * ((0.4 / 2) * (1 + (tanh( ((((( (log10( ec_Ca)) - -0.4) * (ec_Vm - -80.8)) -53.3) / ((1.32e-3 * (ec_Vm + (53.3*( (log10(ec_Ca))--0.4))--80.8)^2 )+0.3)))))) + (0.6 / 2) * (1 + (tanh((((log10( ec_Ca)) - -0.28) / 0.389))))) + 955 * (ec_Vm - -31.1))) - Vm_ht*(ec_Vm - smc_Vm);                                                                 
dxdt(9) = JPLC - 0.1 * ec_IP3 - IP3_ht*(ec_IP3 - smc_IP3);                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                