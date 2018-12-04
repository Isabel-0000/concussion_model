% Isabel Jellen
% BMED 430
% Term Project
% Concussion Indication Model Validation
% File: curvefit.m

% Parameter: linear and rotational acceleration arrays
% Returns: array of experimentally derived concussion probablities
function experimental_cp = curvefit(data)
   
    experimental_cp = zeros(length(data.linear), 1);
    for i = 1:length(data.linear)
        a = data.linear(i);
        alpha = data.rotational(i);
        num = -(-10.2+0.0433.*a+0.000873.*alpha-0.00000092.*a.*alpha);
        experimental_cp(i) = 1/(1+exp(num));
    end
end