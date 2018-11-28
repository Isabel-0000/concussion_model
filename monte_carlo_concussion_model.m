% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: monte_carlo_concussion_model.m

data = generate_data();
concussion_probability_curvefit = curvefit(data)
concussion_probability_model = model(data)
%plot(data, concussion_probability_curvefit, data, concussion_probability_model);