% Define parameters
k_prod = 1;                % ABA production rate constant (concentration/time)
k_deg = 0.1;               % ABA degradation rate constant (1/time)
k_on = 1e4;                % ABA binding to receptor, association rate constant (M^-1 s^-1)
k_off = 1e-3;              % ABA dissociating from receptor, dissociation rate constant (s^-1)
k_half = 1e-2;             % Half-maximal concentration for the Hill equation
n = 2;                     % Cooperativity coefficient in ligand binding for Hill equation
k_trans = 1;               % Transcription rate constant for sensor protein mRNA (concentration/time)
k_trans_deg = 0.1;         % Degradation rate constant for sensor protein mRNA (1/time)
k_trans_syn = 1;           % Translation rate constant for sensor protein synthesis (concentration/time)
k_deg_prot = 0.1;          % Degradation rate constant for sensor protein (1/time)
k_response = 0.25;          % Rate constant for bioluminescence response (1/time)
R_total = 1;               % Total receptor concentration (concentration)

% Define initial conditions
C_ABA_init = 1.5:5:50;       % Initial ABA concentration range (0 to 50 µM)
k_prod_values = linspace(0.5, 1.5, 11);
tspan = [0 100];
threshold = 0.5;
C_ABA_R_init = 0;          % Initial concentration of ABA-receptor complex
mRNA_init = 0.1;           % Small initial concentration of sensor protein mRNA
Protein_init = 0.1;        % Small initial concentration of sensor protein
I_Bioluminescence_init = 0; % Initial bioluminescence response

% Initialize storage for results
final_bioluminescence = zeros(length(C_ABA_init), length(k_prod_values));
detection_time = zeros(length(C_ABA_init), length(k_prod_values));

% Solve the ODEs for different initial ABA concentrations and k_prod values
for i = 1:length(C_ABA_init)
    for j = 1:length(k_prod_values)
        % Current initial conditions and k_prod value
        C_ABA = C_ABA_init(i);
        k_prod = k_prod_values(j);
        initial_conditions = [C_ABA, 0, mRNA_init, Protein_init, I_Bioluminescence_init];  % Updated initial conditions
        
        % Parameter struct with stochasticity added
        params = struct('k_prod', k_prod + randn * 0.05, 'k_deg', k_deg + randn * 0.01, ...
                        'k_on', k_on, 'k_off', k_off, 'k_half', k_half, 'n', n, ...
                        'k_trans', k_trans, 'k_trans_deg', k_trans_deg, ...
                        'k_trans_syn', k_trans_syn, 'k_deg_prot', k_deg_prot, ...
                        'k_response', k_response + randn * 0.05, 'R_total', R_total);
        
        % Solve ODEs
        [T, Y] = ode15s(@(t, Y) aba_biosensor_odes(t, Y, params), tspan, initial_conditions);
        
        % Store final bioluminescence response
        final_bioluminescence(i, j) = Y(end, 5);
        
        % Calculate detection time
        idx = find(Y(:, 5) >= threshold, 1);
        if ~isempty(idx)
            detection_time(i, j) = T(idx);
        else
            detection_time(i, j) = NaN;  % If threshold is not reached within the time span
        end
    end
end

% Create 3D plot for final bioluminescence response
[C_ABA_grid, k_prod_grid] = meshgrid(C_ABA_init, k_prod_values);
figure;
surf(C_ABA_grid, k_prod_grid, final_bioluminescence');
xlabel('Initial ABA Concentration (\muM)');
ylabel('ABA Production Rate (k_{prod})');
zlabel('Final Bioluminescence Response');
title('3D Plot of Bioluminescence Response');
colorbar;

% Define the system of ODEs
function dYdt = aba_biosensor_odes(t, Y, params)
    % Unpack parameters
    k_prod = params.k_prod;
    k_deg = params.k_deg;
    k_on = params.k_on;
    k_off = params.k_off;
    k_half = params.k_half;
    n = params.n;
    k_trans = params.k_trans;
    k_trans_deg = params.k_trans_deg;
    k_trans_syn = params.k_trans_syn;
    k_deg_prot = params.k_deg_prot;
    k_response = params.k_response;
    R_total = params.R_total;
    
    % Unpack state variables
    C_ABA = Y(1);
    C_ABA_R = Y(2);
    mRNA = Y(3);
    Protein = Y(4);
    I_Bioluminescence = Y(5);
    
    % Define ODEs
    dC_ABA_dt = k_prod - k_deg * C_ABA - k_on * C_ABA * Protein + k_off * C_ABA_R;
    dC_ABA_R_dt = k_on * C_ABA * Protein - k_off * C_ABA_R;
    dmRNA_dt = k_trans * (C_ABA^n / (k_half^n + C_ABA^n)) - k_trans_deg * mRNA;
    dProtein_dt = k_trans_syn * mRNA - k_deg_prot * Protein;
    dI_Bioluminescence_dt = k_response * (C_ABA_R^n / (k_half^n + C_ABA_R^n));
    
    % Pack derivatives
    dYdt = [dC_ABA_dt; dC_ABA_R_dt; dmRNA_dt; dProtein_dt; dI_Bioluminescence_dt];
end
