%Function to evaluate polynomial gradient as sum of basis functions
%Total result is lapu = SUM(i, c(i)*LapPhi(i,x))

%NOTE: Matlab indexing starts from 1
function lapu = EvalPolySecondDerivative(c, x)
    lapu = 0;
    for i = 1:size(c,1)
        lapu = lapu + c(i,1)*SecondDerivativePolyBasis1D(i,x);
    end
end

