% Isabel Jellen
% Runner for concussion modeling project
% Validation, propagation of error, and timing 


clear all;

d = generate_data(1000);

tic;
cp2 = model(d);
toc
tic;
cp = curvefit(d);
toc


model_time = zeros(30, 1);

for i = 1: 30
   d = generate_data(i);
   tstart = tic;
   model(d);
   model_time(i) = toc(tstart);   
end
datasize = linspace(1, 30, 30);
[fit_linear, S] = polyfit(datasize.',model_time,1);
p = polyval(fit_linear, datasize.')
figure(6)
scatter(datasize, model_time);
hold on;
plot(p);
title("Model Time based on Dataset Size");
ylabel("Model time (s)");
xlabel("Dataset size");

hold off;



t_final = d.duration(1) + 85*10^-3;
time_vector = linspace(0, t_final, 100);

shape_fxn = d.alpha_of_t_x(1,:)./d.rotational_x(1);
figure(4)
hold on;
plot(time_vector, shape_fxn);
title("Head Acceleration Shape Function");
ylabel("Peak acceleration scalar multiple");
xlabel("Time after impact (s)");
hold off;

    

[h,p,ci,stats] = ttest2(cp, cp2)
rms = sqrt(mean((cp - cp2).^2))
rms_max = sqrt(max(abs(cp-cp2))^2);

figure(2)

scatter(d.rotational, cp);
hold on;
title("Concussion Probability based on Rotational Acceleration");
ylabel("Concussion Probability");
xlabel("Rotational Acceleration Magnitude (rad/s^2)");
scatter(d.rotational, cp2);
legend({'Curve Fit', 'Model'},'Location', 'southeast');
hold off;
alpha(0.5);

figure(3)
histogram(d.rotational);
hold on;
xlabel("Rotational Acceleration Magnitude (rad/s^2)");
ylabel("Number of Impacts");
title("Histogram of Rotational Acceleration Data");
hold off;


figure(5)
scatter(d.rotational, cp-cp2);
title("Residual in Concussion Probability based on Rotational Acceleration");
ylabel("Residual between curve fit and model");
xlabel("Rotational Acceleration Magnitude (rad/s^2)");

% Change rotational acceleration by 5% "error"
noisy_data = d;
random_noise = (1.05 - 0.95).*rand(1000, 1) + 0.95;
noisy_data.linear = noisy_data.linear.*random_noise;
noisy_data.rotational = noisy_data.rotational.*random_noise;
noisy_data.alpha_of_t_x = noisy_data.alpha_of_t_x.*random_noise;
noisy_data.alpha_of_t_y = noisy_data.alpha_of_t_y.*random_noise;
noisy_data.omega_p = noisy_data.omega_p.*random_noise;
noisy_data.omega_p_x = noisy_data.omega_p_x.*random_noise;
noisy_data.omega_p_y = noisy_data.omega_p_y.*random_noise;

cp2_bad_data = model(noisy_data);
cp_bad_data = curvefit(noisy_data);

rms_c = sqrt(mean((cp - cp_bad_data).^2))
rms_m = sqrt(mean((cp2 - cp2_bad_data).^2))
[h_5c,p_5c,ci_5c,stats_5c] = ttest2(cp, cp_bad_data)
[h_5m,p_5m,ci_5m,stats_5m] = ttest2(cp2, cp2_bad_data)

figure(4)
scatter(d.rotational, cp);
hold on;
title("Concussion Probability based on Rotational Acceleration");
ylabel("Concussion Probability");
xlabel("Rotational Acceleration Magnitude (rad/s^2)");
scatter(d.rotational, cp2);
scatter(d.rotational, cp_bad_data);
scatter(d.rotational, cp2_bad_data);
legend({'Curve Fit', 'Model', 'Curve Fit - 5% error', 'Model - 5% error'},'Location', 'southeast');
hold off;
alpha(0.5);




