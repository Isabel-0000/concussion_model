% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: generate_data.m

% "For eight players there were 347 total impacts with an average of 21.5 ±
%       19.7 g, 11.5 ± 33.4 HIC, 16.7 ± 51.2 GSI, 769.9 ± 1082.7 rad/s2 rotation about the x-axis and
%       1382.8 ± 1547.3 rad/s2 rotation about the y-axis" [1]
% [1] https://www-nrd.nhtsa.dot.gov/pdf/bio/proceedings/2003_31/31-6.pdf


% Create normal distributions of 1000 samples of 
%       impacts generated from experimental data
% Returns: linear and rotational acceleration arrays
function data = generate_data()
    g = 9.81; % [m/s^2]
    % Linear acceleration
    pdfs(1) = makedist('Normal', 'mu', 21.5 * g, 'sigma', 19.7 * g); % [m/s^2] 
    % Dominant rotation about the y-axis
    pdfs(2) = makedist('Normal', 'mu', 1382.8, 'sigma', 1547.3); % [rad/s^2]
    % Duration of impact: chosen based upon wayne state tolerance curve
    % seen here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2790598/
    pdfs(3) = makedist('Normal', 'mu', 6*10^-3, 'sigma', 2*10^-3);
    data.linear = abs(random(pdfs(1), 1, 1));
    data.rotational = abs(random(pdfs(2), 1, 1));
    data.duration = abs(random(pdfs(3), 1, 1));
    
    data.alpha_of_t = head_acceleration_shape_function(data.rotational, data.duration); % time varying rotational acceleration
    data.omega_p = peak_angular_velocity(data.rotational, data.duration);
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
        
        %{
        syms t;
        shape_fxn = piecewise (0 <= t < d/2, (a/(d/2))*t, ...
                        d/2 <= t < d, 2*a-(a*2/(d))*t, ...
                        d <= t < d + 35/2*10^-3, d*a*.4296/(35/2*10^-3)-((a*.4286)/(35/2*10^-3))*t, ...
                        d + 35/2*10^-3 <= t < d + 35*10^-3, ((a*.4286)/(35/2*10^-3))*t - 2*0.4296*a - d*0.4296*a/(35/2*10^-3), ...
                        d + 35*10^-3 <= t <= d + 85*10^-3, 0);
        %}
    end
end    

function omega_peak = peak_angular_velocity(peak_alpha, duration)
    % Integral of the triangle seen here
    % https://academic.oup.com/milmed/article/183/suppl_1/339/4959982
    omega_peak = peak_alpha.*duration;  
end