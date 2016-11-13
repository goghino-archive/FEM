function U0 = u0(x,y,z)
% u0 is Dirichlet's boundary function of the poisson's equation
% and also the exact analytical solution

U0 = x .* exp(- (y-1).^2 .* (z-1).^2 );

end