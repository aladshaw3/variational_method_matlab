%Function for the simple 1-D polynomial basis
%Basis function is phi = x^(i-1)

%NOTE: Matlab indexing starts from 1
function phi = PolyBasis1D(i, x)
    if i <= 1
        phi = 1;
    else
        phi = x^(i-1);
    end
end