% Plot the "bifurcation diagram" of JPLC vs smc_Ca for the 3 cases.
% Case 1: Vm_ht = 50, Ca_ht = 0, IP3_ht = 0.05
% Case 2/3: Vm_ht = 50, Ca_ht = 0.05, IP3_ht = 0.05
% Case 4: Vm_ht = 0, Ca_ht = 0, IP3_ht = 0.05

% Just press run!

clear;
odeopts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.1, 'Vectorized', 1);

% Case 1
% Oscillations between JPLC = [0.3668, 0.9046]
Vm_ht = 50;
Ca_ht = 0;
IP3_ht = 0.05;

JPLC_range = 0:0.005:1;
i = 1;
max_val_1 = zeros(size(JPLC_range));
min_val_1 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) CoupledCellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht);
    [t, x] = ode15s(f, 0:0.01:1000, [0.1,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca = x((length(x)*0.9):length(x),1); 
    
    max_val_1(i) = max(smc_Ca);
    min_val_1(i) = min(smc_Ca);
    i = i+1
end

figure(20);
plot(JPLC_range, max_val_1, 'b', JPLC_range, min_val_1, 'b');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 1: Vm_{ht} = 50, Ca_{ht} = 0, IP3_{ht} = 0.05');

% Case 2/3
% Oscillations between JPLC = [0.2688, 0.5781]
Vm_ht = 50;
Ca_ht = 0.05;
IP3_ht = 0.05;

j = 1;
max_val_2 = zeros(size(JPLC_range));
min_val_2 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) CoupledCellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht);
    [t, x] = ode15s(f, 0:0.01:1000, [0.1,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca = x((length(x)*0.9):length(x),1); 
    
    max_val_2(j) = max(smc_Ca);
    min_val_2(j) = min(smc_Ca);
    j = j+1
end

figure(3);
plot(JPLC_range, max_val_2, 'r', JPLC_range, min_val_2, 'r');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 2/3: Vm_{ht} = 50, Ca_{ht} = 0.05, IP3_{ht} = 0.05');

% Case 4
% Oscillations between JPLC = [0.2939, 0.6534]
Vm_ht = 0;
Ca_ht = 0;
IP3_ht = 0.05;

k = 1;
max_val_4 = zeros(size(JPLC_range));
min_val_4 = zeros(size(JPLC_range));

for JPLC = JPLC_range
    f = @(t,x) CoupledCellsODE(t,x,JPLC,Vm_ht,Ca_ht,IP3_ht);
    [t, x] = ode15s(f, 0:0.01:1000, [0.1,0.1,-60,0.1,0.1,0.12,0.12,-62,0.12], odeopts);

    % Remove the first 90% of solution (transient behaviour)
    smc_Ca = x((length(x)*0.9):length(x),1); 
    
    max_val_4(k) = max(smc_Ca);
    min_val_4(k) = min(smc_Ca);
    k = k+1
end

figure(4);
plot(JPLC_range, max_val_4, 'g', JPLC_range, min_val_4, 'g');
xlabel('J_{PLC}');
ylabel('smc_{Ca}');
title('Case 4: Vm_{ht} = 0, Ca_{ht} = 0, IP3_{ht} = 0.05');

% Save data in file in case of power outage etc.
save('CoupledCellsODEbifurcation.mat');
