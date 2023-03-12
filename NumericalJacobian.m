function J = NumericalJacobian(fun,x)
    tol = sqrt(eps);    % Use Sqrt of machine precision as a tolerance
    F = fun(x);
    N = length(x);
    M = length(F);
    J = zeros(M,N);
    for i=1:N
        dx = x;
        dx(i) =  x(i) + tol;
        J(:, i) = (fun(dx) - F)/tol;
    end
end

