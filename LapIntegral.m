%Function for the Laplacian Integral of the simple 1-D polynomial basis
%   Basis function is phi = x^(i-1)
%   Second derivative is lapphi = (j-1)*(j-2)*x^(j-3)
%   
%   Function will perform integration of phi*lapphi
%   with lower and upper limits

%NOTE: Matlab indexing starts from 1
function Aij = LapIntegral(i, j, lower, upper)
    if j <= 2
        Aij = 0;
    else
        exp = (i-1)+(j-1)-1;
        pre = (j-1)*(j-2);
        Aij = pre * ((upper^exp /exp) - (lower^exp/exp));
    end
end

