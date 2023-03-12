function Jv = JacobianOperator(fun,v,x)
    tol = sqrt(eps);    % Use Sqrt of machine precision as a tolerance
    Jv= (fun(x+tol*v) - fun(x))/tol;
end

