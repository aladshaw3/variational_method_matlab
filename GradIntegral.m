%Function for the Gradient Integral of the simple 1-D polynomial basis
%   Basis function is phi = x^(i-1)
%   First derivative is gradphi = (j-1)*x^(j-2)
%   
%   Function will perform integration of phi*gradphi
%   with lower and upper limits

%NOTE: Matlab indexing starts from 1
function Bij = GradIntegral(i, j, lower, upper)
    if j <= 1
        Bij = 0;
    else
        exp = (i-1)+(j-1);
        pre = (j-1);
        Bij = pre * ((upper^exp /exp) - (lower^exp/exp));
    end
end

