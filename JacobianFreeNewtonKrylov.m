function [x,nl_iter] = JacobianFreeNewtonKrylov(fun,x0)
    tol = sqrt(eps);    % Use Sqrt of machine precision as a tolerance
    x = x0;
    F = fun(x0);
    nl_iter=0;

    while (norm(F) > tol)
        Jop = @(F) JacobianOperator(fun,F,x);
        s = gmres(Jop,-F);
        x = x + s;
        F = fun(x);
        nl_iter=nl_iter+1;
    end
end

