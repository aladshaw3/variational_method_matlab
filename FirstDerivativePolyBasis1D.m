%Function for the 1st derivative of the simple 1-D polynomial basis
%   Basis function is phi = x^(i-1)
%   First derivative is gradphi = (i-1)*x^(i-2)

%NOTE: Matlab indexing starts from 1
function gradphi = FirstDerivativePolyBasis1D(i, x)
    if i <= 1
        gradphi = 0;
    else
        gradphi = (i-1)*x^(i-2);
    end
end

