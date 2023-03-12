%Function for the Overlap Integral of the simple 1-D polynomial basis
%   Basis function is phi = x^(i-1)
%   
%   Function will perform integration of phi*phi
%   with lower and upper limits

%NOTE: Matlab indexing starts from 1
function Oij = OverlapIntegral(i, j, lower, upper)
    exp = (i-1)+(j-1)+1;
    Oij = ((upper^exp /exp) - (lower^exp/exp));
end

