% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: curvefit.m

% Parameter: linear and rotational acceleration arrays
% Returns: array of experimentally derived concussion probablities
function experimental_cp = curvefit(data)
    a = data.linear;
    alpha = data.rotational;
    num = -(-10.2+0.0433.*a+0.000873.*alpha-0.00000092.*a.*alpha);
    experimental_cp = 1/(1+exp(num));
end