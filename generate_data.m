% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: generate_data.m

% "For eight players there were 347 total impacts with an average of 21.5 �
%       19.7 g, 11.5 � 33.4 HIC, 16.7 � 51.2 GSI, 769.9 � 1082.7 rad/s2 rotation about the x-axis and
%       1382.8 � 1547.3 rad/s2 rotation about the y-axis" [1]
% [1] https://www-nrd.nhtsa.dot.gov/pdf/bio/proceedings/2003_31/31-6.pdf


% Create normal distributions of 1000 samples of 
%       impacts generated from experimental data
% Returns: linear and rotational acceleration arrays
function data = generate_data(n)
    g = 9.81; % [m/s^2]
    lin_mean = 21.5;
    lin_sd = 19.7;
    x_ang_mean = 769.9;
    y_ang_mean = 1382.2;
    %x_ang_mean = 561;
    %y_ang_mean = 625;
    z_ang_mean = 0; % 1106
   
    % Linear acceleration
    pdfs(1) = makedist('Normal', 'mu', lin_mean * g, 'sigma', lin_sd * g); % [m/s^2] 
   
    % Duration of impact: chosen based upon wayne state tolerance curve
    % seen here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2790598/
    %   - Chose a SD of 0 to minimize variability in data in order to get
    %     the cleanest curve
    pdfs(3) = makedist('Normal', 'mu', 6*10^-3, 'sigma', 0*10^-3);
    
    % Create data arrays from pdfs
    data.linear = abs(random(pdfs(1), n, 1));
    % For rotational acceleration:
    %   - take pythagorean sum of mean angular accelerations given
    %   - scale linear data by this sum divided by the mean linear
    %     acceleration
    data.rotational = data.linear.*sqrt(x_ang_mean^2 + y_ang_mean^2 + z_ang_mean^2)/lin_mean;
    data.duration = abs(random(pdfs(3), n, 1));
    data.omega_p = peak_angular_velocity(data.rotational, data.duration);
    data.rotational_x = data.linear.*x_ang_mean/lin_mean;
    data.rotational_y = data.linear.*y_ang_mean/lin_mean;
    data.rotational_z = data.linear.*z_ang_mean/lin_mean;
    
    % time varying rotational acceleration
    data.alpha_of_t_x = head_acceleration_shape_function(data.rotational_x, data.duration); 
    data.alpha_of_t_y = head_acceleration_shape_function(data.rotational_y, data.duration);
    data.alpha_of_t_z = head_acceleration_shape_function(data.rotational_z, data.duration);
    data.omega_p_x = peak_angular_velocity(data.rotational_x, data.duration);
    data.omega_p_y = peak_angular_velocity(data.rotational_y, data.duration);
    data.omega_p_z = peak_angular_velocity(data.rotational_z, data.duration);
end

% Create a head acceleration shape function based upon the shape given here:
% https://academic.oup.com/milmed/article/183/suppl_1/339/4959982
function shape_fxn = head_acceleration_shape_function(peak_alpha, duration)
    for j = 1: length(peak_alpha)
        d = duration(j);
        a = peak_alpha(j);
        
        t = linspace(0, d + 85*10^-3, 100);
        for i = 1: length(t)
            if(t(i) < d/2)
                y(i) = (a/(d/2))*t(i);
            elseif(t(i) < d)
                y(i) = 2*a-(a*2/(d))*t(i);
            elseif(t(i) < d + 35/2*10^-3)
                y(i) = d*a*.4296/(35/2*10^-3)-((a*.4286)/(35/2*10^-3))*t(i);
            elseif(t(i) < d + 35*10^-3)
                y(i) = ((a*.4286)/(35/2*10^-3))*t(i) - (d + 35*10^-3)*((a*.4286)/(35/2*10^-3));
            else
                y(i) = 0;
            end
            shape_fxn(j,i) = y(i);
        end
        
    end
end    

function omega_peak = peak_angular_velocity(peak_alpha, duration)
    % Integral of the triangle seen here
    % https://academic.oup.com/milmed/article/183/suppl_1/339/4959982
    omega_peak = peak_alpha.*duration;  
end