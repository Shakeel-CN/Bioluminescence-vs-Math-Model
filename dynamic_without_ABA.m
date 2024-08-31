% Define symbolic variables
syms SNRK2(t) PP2C(t) MAPK(t)

% Interaction constants
k_interact_SNRK2 = 0.1;
k_interact_PP2C = 0.2;
k_interact_MAPK = 0.3;

% Set initial conditions
SNRK2_0 = 1.0; % Initial SNRK2 concentration
PP2C_0 = 0.5;  % Initial PP2C concentration
MAPK_0 = 0.2;  % Initial MAPK concentration

% Time span for simulation
tspan = 0:1:50;

% Define the ODE system for protein interactions
odeSystem = @(t, y) [
    k_interact_SNRK2 * y(2) - k_interact_PP2C * y(3);  % Equation for dPP2C/dt
    -k_interact_SNRK2 * y(1) + k_interact_MAPK * y(3);  % Equation for dMAPK/dt
    -k_interact_MAPK * y(2);  % Equation for dSNRK2/dt
];

% Solve the ODE numerically
[t, y] = ode45(odeSystem, tspan, [SNRK2_0, PP2C_0, MAPK_0]);

% Plotting SNRK2, PP2C, and MAPK interactions over time
figure;
plot(t, y(:, 1), '-');  % SNRK2
hold on;
plot(t, y(:, 2), '--'); % PP2C
plot(t, y(:, 3), '-.'); % MAPK

xlabel('Time (sec)');
ylabel('Concentration (μmol)');
legend('SNRK2', 'PP2C', 'MAPK');
title('Interaction of Proteins without ABA Over Time');
hold off;
% Plotting SNRK2, PP2C, and MAPK interactions over time
figure('Units', 'normalized', 'OuterPosition', [0 0 0.8 0.8]);
plot(t, y(:, 1), '-');  % SNRK2
hold on;
plot(t, y(:, 2), '--'); % PP2C
plot(t, y(:, 3), '-.'); % MAPK
xlabel('Time (sec)');
ylabel('Concentration (μmol)');
legend('SNRK2', 'PP2C', 'MAPK');
title('Interaction of Proteins without ABA Over Time');
hold off;