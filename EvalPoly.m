%Function to evaluate polynomial as sum of basis functions
%Total result is u = SUM(i, c(i)*Phi(i,x))

%NOTE: Matlab indexing starts from 1
function u = EvalPoly(c, x)
    u = 0;
    for i = 1:size(c,1)
        u = u + c(i,1)*PolyBasis1D(i,x);
    end
end

