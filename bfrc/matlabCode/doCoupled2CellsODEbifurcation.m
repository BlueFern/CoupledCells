% Plot the "bifurcation diagram" of JPLC vs smc_Ca_1 and smc_Ca_2 for the 4 cases.
% Dotted line - second smc
% Normal line - first smc
% Both cells should be the same since they have the same JPLC values -
% dotted line is a check

% Just press run!

clear;
odeopts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.1, 'Vectorized', 1);
JPLC_range = 0:0.005:1;


% Case 1
Vm_ht = 50;
Ca_ht = 0;
IP3_ht = 0.05;
Vm_hm_smc = 1000;
Ca_hm_smc = 0.05;
IP3_hm_smc = 0.05;
Vm_hm_ec = 1000;
Ca_hm_ec = 0.05;
IP3_hm_ec = 0;

i = 1;
max_val_1_1 = zeros(size(JPLC_range));
min_val_1_1 = zeros(size(JPLC_range));
max_val_1_2 = zeros(size(JPLC_range));
min_val_1_2 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) Coupled2CellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht,Vm_hm_smc,Ca_hm_smc,IP3_hm_smc,Vm_hm_ec,Ca_hm_ec,IP3_hm_ec);
    [t x] = ode15s(f, 0:0.01:1000, [0.3,0.2,-60,0.1,0.4,0.12,0.12,-62,0.12,0.4,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca_1 = x((length(x)*0.9):length(x),1); 
    smc_Ca_2 = x((length(x)*0.9):length(x),10); 
    
    max_val_1_1(i) = max(smc_Ca_1);
    min_val_1_1(i) = min(smc_Ca_1);
    max_val_1_2(i) = max(smc_Ca_2);
    min_val_1_2(i) = min(smc_Ca_2);
    i = i+1
end

figure(10);
plot(JPLC_range, max_val_1_1, 'b', JPLC_range, min_val_1_1, 'b', JPLC_range, max_val_1_2, ':b', JPLC_range, min_val_1_2, ':b');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 1');



% Case 2
Vm_ht = 50;
Ca_ht = 0.05;
IP3_ht = 0.05;
Vm_hm_smc = 1000;
Ca_hm_smc = 0.05;
IP3_hm_smc = 0.05;
Vm_hm_ec = 1000;
Ca_hm_ec = 0.05;
IP3_hm_ec = 0;

j = 1;
max_val_2_1 = zeros(size(JPLC_range));
min_val_2_1 = zeros(size(JPLC_range));
max_val_2_2 = zeros(size(JPLC_range));
min_val_2_2 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) Coupled2CellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht,Vm_hm_smc,Ca_hm_smc,IP3_hm_smc,Vm_hm_ec,Ca_hm_ec,IP3_hm_ec);
    [t x] = ode15s(f, 0:0.01:1000, [0.3,0.2,-60,0.1,0.4,0.12,0.12,-62,0.12,0.4,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca_1 = x((length(x)*0.9):length(x),1); 
    smc_Ca_2 = x((length(x)*0.9):length(x),10); 
    
    max_val_2_1(j) = max(smc_Ca_1);
    min_val_2_1(j) = min(smc_Ca_1);
    max_val_2_2(j) = max(smc_Ca_2);
    min_val_2_2(j) = min(smc_Ca_2);
    j = j+1
end

figure(20);
plot(JPLC_range, max_val_2_1, 'r', JPLC_range, min_val_2_1, 'r', JPLC_range, max_val_2_2, ':r', JPLC_range, min_val_2_2, ':r');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 2');




% Case 3
Vm_ht = 50;
Ca_ht = 0.05;
IP3_ht = 0.05;
Vm_hm_smc = 1000;
Ca_hm_smc = 0.05;
IP3_hm_smc = 0.05;
Vm_hm_ec = 1000;
Ca_hm_ec = 0.05;
IP3_hm_ec = 0.05;

k = 1;
max_val_3_1 = zeros(size(JPLC_range));
min_val_3_1 = zeros(size(JPLC_range));
max_val_3_2 = zeros(size(JPLC_range));
min_val_3_2 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) Coupled2CellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht,Vm_hm_smc,Ca_hm_smc,IP3_hm_smc,Vm_hm_ec,Ca_hm_ec,IP3_hm_ec);
    [t x] = ode15s(f, 0:0.01:1000, [0.3,0.2,-60,0.1,0.4,0.12,0.12,-62,0.12,0.4,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca_1 = x((length(x)*0.9):length(x),1); 
    smc_Ca_2 = x((length(x)*0.9):length(x),10); 
    
    max_val_3_1(k) = max(smc_Ca_1);
    min_val_3_1(k) = min(smc_Ca_1);
    max_val_3_2(k) = max(smc_Ca_2);
    min_val_3_2(k) = min(smc_Ca_2);
    k = k+1
end

figure(30);
plot(JPLC_range, max_val_3_1, 'k', JPLC_range, min_val_3_1, 'k', JPLC_range, max_val_3_2, ':k', JPLC_range, min_val_3_2, ':k');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 3');






% Case 4
Vm_ht = 0;
Ca_ht = 0;
IP3_ht = 0.05;
Vm_hm_smc = 1000;
Ca_hm_smc = 0.05;
IP3_hm_smc = 0.05;
Vm_hm_ec = 0;
Ca_hm_ec = 0;
IP3_hm_ec = 0.05;

l = 1;
max_val_4_1 = zeros(size(JPLC_range));
min_val_4_1 = zeros(size(JPLC_range));
max_val_4_2 = zeros(size(JPLC_range));
min_val_4_2 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) Coupled2CellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht,Vm_hm_smc,Ca_hm_smc,IP3_hm_smc,Vm_hm_ec,Ca_hm_ec,IP3_hm_ec);
    [t x] = ode15s(f, 0:0.01:1000, [0.3,0.2,-60,0.1,0.4,0.12,0.12,-62,0.12,0.4,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca_1 = x((length(x)*0.9):length(x),1); 
    smc_Ca_2 = x((length(x)*0.9):length(x),10); 
    
    max_val_4_1(l) = max(smc_Ca_1);
    min_val_4_1(l) = min(smc_Ca_1);
    max_val_4_2(l) = max(smc_Ca_2);
    min_val_4_2(l) = min(smc_Ca_2);
    l = l+1
end

figure(40);
plot(JPLC_range, max_val_4_1, 'g', JPLC_range, min_val_4_1, 'g', JPLC_range, max_val_4_2, ':g', JPLC_range, min_val_4_2, ':g');

xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 4');

% Save data in file in case of power outage or matlab crashing etc.
save('Coupled2CellsODEbifurcation.mat');