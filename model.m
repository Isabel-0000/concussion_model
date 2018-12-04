% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: model.m
% Sources: https://academic.oup.com/milmed/article/183/suppl_1/339/4959982
%          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471336/


% Parameter: linear and rotational acceleration arrays
% Returns: array of modeled concussion probablities
function model_cp = model(data)
    % Assume rotation about y axis dominates (see generate_data)
    
    % Step 1: simplified human head FEA
    model_cp = zeros(length(data.omega_p), 1);
    for i = 1:length(data.omega_p)
        alpha_x = data.alpha_of_t_x(i,:);
        alpha_y = data.alpha_of_t_y(i,:);
        omega_p_x = data.omega_p_x(i);
        omega_p_y = data.omega_p_y(i);
        duration = data.duration(i);
        
        params.a = 3.3;
        params.b = 250;
        params.c = -2.2;
        params.d = 74800;
        strain_x = simplified_hh_FEM(alpha_x, params, omega_p_x, duration);
        params.a = -3.0;
        params.b = -230;
        params.c = 3.6;
        params.d = -67320;
        strain_y = simplified_hh_FEM(alpha_y, params, omega_p_y, duration);
        
        params.a = 2.7;
        params.b = 150;
        params.c = 0;
        params.d = 59840;
        strain_z = simplified_hh_FEM(alpha_y, params, omega_p_y, duration);
        
        nr_strain_x = strain_node_ranvier(strain_x, duration);
        nr_strain_y  = strain_node_ranvier(strain_y, duration);
        nr_strain_z  = strain_node_ranvier(strain_z, duration);
        nr_strain = sqrt(nr_strain_x^2 + nr_strain_y^2 + nr_strain_z^2);
        delta_action_potential = axon_signaling(nr_strain);
        model_cp(i) = dose_response(delta_action_potential);
    end
    
end

% Model 1: simplification of FEM of a head to find strain
function strain = simplified_hh_FEM(alpha, params, omega_p, duration)
    
    t_final = duration + 85*10^-3;
    time_vector = linspace(0, t_final, 100);
    tspan = [0 t_final];
    y0 = [0 0];
    [t,y] = ode45(@(t, y) second_order_ode(t, y, alpha, omega_p, params, time_vector), tspan, y0);
    strain = y(:,1);
end

function dydt = second_order_ode(t, y, alpha, omega_p, params, time_vector)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    f = interp1(time_vector, alpha, t);
    dydt(2) = (1/params.a)*...
        (-f - ...
        (params.b + ...
        params.c*omega_p)*...
        y(2) - ...
        params.d*y(1));
end

% Model 2: translate strain to strain at node of ranvier (this was super
% fun - like 9 equations and 9 unknowns that were simplified into 3 odes
% kind of fun)
function node_of_ranvier_strain = strain_node_ranvier(input_strain, duration)
    % The following constants are from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471336/ 
    L_node = 1*10^-6; % 1 micron
    L_internode = 125*10^-6; % 50-200 microns
    % From "Electrical Properties of the axon"
    A_node = 7.85*10^-11;
    A_internode = 7.85*10^-11;
    A_myelin = 7.54*10^-11;
    
    E1 = 19.9*10^3;
    E2 = 0.42*10^3;
    nu1 = 2.256*10^6;
    E3 = 50*10^3;
    nu3 = 1*10^3;
    
    % k means spring constant, d means damping constant
    % See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5471336/ for
    % viscoelastic circuit model
    
    c.k11 = E1*A_internode/L_internode; 
    c.k12 = E2*A_internode/L_internode;
    c.d1 = nu1*A_internode/L_internode;
    c.k21 = E1*A_node/L_node;
    c.k22 = E2*A_node/L_node;
    c.d2 = nu1*A_node/L_node;
    c.k3 = E3*A_myelin/L_internode;
    c.d3 = nu3*A_myelin/L_internode;
    
    t_final = duration + 85*10^-3;
    tspan = [0 t_final];
    y0 = [0 0 0]; % All zero ICs
    % Get strains at each damper e_n1, e_n2, e_n3
    [t,y] = ode45(@(t, y) odefun_micromech(t, y, c), tspan, y0);

    % Back out strains at each individual spring
    % Except e_k3 which we will figure out doesn't really matter
    e_k22 = y(2:end,2);
    e_k12 = y(2:end,1);
    e_k11 = c.d3*diff(y(:,3))+ y(3) - y(1);
    e_k21 = (1/c.k21)*(y(1)-y(3) + (1 + c.k11)*e_k11);
    %e_k21 = (c.k22.*y(2:end,2) + c.d3.*diff(y(:,2)))./c.k21;
    %e_k11 = (c.k21.*e_k21 + c.k3.*y(2:end,1)./y(2:end,2))./(c.k11 - c.k3./y(2:end,2));
    
    % Plug individual strains into final eqn to get final strain
    % First, trim input
    time = linspace(0, t_final, length(y(2:end,2)));
    input_time = linspace(0, t_final, length(input_strain));
    input_strain = interp1(input_time,input_strain,time).';
    % Second, find strain
    e_out = input_strain - e_k11 - e_k12 - e_k21 - e_k22;
    % Third, find max and return it
    node_of_ranvier_strain = max(abs(e_out));
end

function dydt = odefun_micromech(t,y,c)
    dydt = zeros(3,1);
    dydt(1) = (1/c.d1)*(c.k22*y(3) + c.d2*dydt(2) - c.d3*dydt(3) - c.k12*y(1));
    dydt(2) = (1/c.d2)*(y(1) - y(3) + (1 + c.k11)*(c.d3*dydt(3) + y(3) - y(1)) - c.k22*y(3));
    dydt(3) = (1/c.d3)*((c.k3*(y(1) - y(2)) + y(3) - y(1))/(1-c.k3) + y(1) - y(3));
end


% Model 3: Axon signaling model
function delta_ap = axon_signaling(strain)
    if(strain > 0.73)
        delta_ap = 1;
    else
        delta_ap = 1.5826*strain^2 + 0.2138*strain;
    end
end

% Model 5: Dose response model
function cp = dose_response(NIM)
    % Determined from curve fit of experimental data in this paper
    % https://academic.oup.com/milmed/article/183/suppl_1/339/4959982
    if(NIM == 1)
        cp = 1;
    else
        w = log(NIM/(1-NIM));
        x = 7.6182 + 2.4587*w;
        cp = exp(x)/(1 + exp(x));
    end
end
    
