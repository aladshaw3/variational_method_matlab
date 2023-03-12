%Main program

%Create basis size and initial guess
N = 5;
x0 = zeros(N+2,1);

%Create object to perform work and set parameters
test = DiffRxn();
test.Length = 1;
test.BasisSize = N;
test.ReactionCoef = 3;
test.DiffusionCoef = 1;
test.BoundaryVal = 1;

% Call solver and store solution in x
%[x, iter] = NewtonMethod(@test.Residual, x0);
[x, iter] = JacobianFreeNewtonKrylov(@test.Residual, x0);
%[x,fval,exitflag,output] = fsolve(@test.Residual, x0);iter=output.iterations;
fprintf(' NL iterations = %i\n', num2str(iter))


%sol = eval('gmres(J,x)')

% Use solution to plot the polynomial curve
[u, z] = test.Evaluate(x);

% Produce exact solution
[u_e, z_e] = test.ExactSoln();

% Plot results
plot(z, u, '*', z_e, u_e)
