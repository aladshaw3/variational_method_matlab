%Function to evaluate polynomial gradient as sum of basis functions
%Total result is gradu = SUM(i, c(i)*GradPhi(i,x))

%NOTE: Matlab indexing starts from 1
function gradu = EvalPolyFirstDerivative(c, x)
    gradu = 0;
    for i = 1:size(c,1)
        gradu = gradu + c(i,1)*FirstDerivativePolyBasis1D(i,x);
    end
end

