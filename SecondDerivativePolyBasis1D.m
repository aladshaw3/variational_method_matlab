%Function for the 2nd derivative of the simple 1-D polynomial basis
%   Basis function is phi = x^(i-1)
%   First derivative is gradphi = (i-1)*x^(i-2)
%   Second derivative is lapphi = (i-1)*(i-2)*x^(i-3)

%NOTE: Matlab indexing starts from 1
function lapphi = SecondDerivativePolyBasis1D(i, x)
    if i <= 2
        lapphi = 0;
    else
        lapphi = (i-1)*(i-2)*x^(i-3);
    end
end

