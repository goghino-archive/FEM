function U0 = u0(X, Y)
% u0 is Dirichlet's boundary function of the poisson's equation
% and also the exact analytical solution

U0 = 10*X + tanh(10*X - 10);

end