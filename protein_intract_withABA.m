% Define symbolic variables
syms A(t) SNRK2(t) PP2C(t) MAPK(t)

% Parameters for the ODE system
K1 = 0.1;   % Rate constant for ABA decay
K2 = 0.5;   % Rate constant for bacterial response
K3 = 0.2;   % Rate constant for bioluminescence
K = 1e-2;   % Half-maximal concentration for the Hill equation
n = 0.5;    % Hill coefficient for the Hill equation

% Interaction constants
k_interact_SNRK2 = 0.1;
k_interact_PP2C = 0.2;
k_interact_MAPK = 0.3;

% Define the Hill equation
hillEquation = A^n / (K^n + A^n);

% Set initial conditions
ABA_0 = 1.0; % Initial ABA concentration
SNRK2_0 = 1.0; % Initial SNRK2 concentration
PP2C_0 = 0.5;  % Initial PP2C concentration
MAPK_0 = 0.2;  % Initial MAPK concentration

% Time span for simulation
tspan = 0:1:50;

% Define the extended ODE system with dynamic concentration changes for SNRK2, PP2C, and MAPK
odeSystem = @(t, y) [
    -K1 * y(1);                                       % Equation for dA/dt
    +K2 * double(subs(hillEquation, A, y(1)));   % Equation for dSNRK2/dt
    k_interact_SNRK2 * y(2) - k_interact_PP2C * y(3);  % Equation for dPP2C/dt
    -k_interact_SNRK2 * y(2) + k_interact_MAPK * y(3);  % Equation for dMAPK/dt
];

% Solve the ODE numerically
[t, y] = ode45(odeSystem, tspan, [ABA_0, SNRK2_0, PP2C_0, MAPK_0]);

% Plotting ABA, SNRK2, PP2C, and MAPK interactions over time
figure;
plot(t, y(:, 1), '-');  % ABA
hold on;
plot(t, y(:, 2), '--'); % SNRK2
plot(t, y(:, 3), '-.'); % PP2C
plot(t, y(:, 4), ':');  % MAPK

xlabel('Time (sec)');
ylabel('Concentration (μmol)');
legend('ABA', 'SNRK2', 'PP2C', 'MAPK');
title('Interaction of Proteins with ABA Over Time');
hold off;