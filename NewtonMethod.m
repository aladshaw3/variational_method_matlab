function [x,nl_iter] = NewtonMethod(fun,x0)
    tol = sqrt(eps);    % Use Sqrt of machine precision as a tolerance
    x = x0;
    F = fun(x0);
    nl_iter=0;

    while (norm(F) > tol)
        J = NumericalJacobian(fun,x);
        % J*s = -F --> s = -J\F
        x = x - J\F;
        F = fun(x);
        nl_iter=nl_iter+1;
    end
end

