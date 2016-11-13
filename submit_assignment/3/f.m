function F = f(X, Y)
% f is RHS of the poisson's equation with vector arguments

%syms x
%f = -(10*x + tanh(10*x - 10))
%diff(f,x,2)

%F = 20.*tanh(10*X - 10).*(10.*tanh(10*X - 10).^2 - 10);
F = -20.*tanh(10*X - 10).*(10*tanh(10*X - 10).^2 - 10);

end