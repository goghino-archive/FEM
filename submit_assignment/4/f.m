function F = f(X, Y, Z)
% f is RHS of the poisson's equation with vector arguments

%% Solution of f given u in Poisson equation:
% - \delta u(x,y,z) = f(x,y,z)
% for Dirichlet and Neumann BC:
% \partial u / \partial n = 0 for y=1, z=1
% u(x,y,z) = u0(x,y,z) elsewhere

% syms x;
% syms y;
% syms z;
% u = @(x,y,z) x * exp(- (y-1)^2 * (z-1)^2 );
% 
% u_xx = diff(u,x,2);
% u_yy = diff(u,y,2);
% u_zz = diff(u,z,2);
% 
% f = -(u_xx + u_yy + u_zz);
% 
% %verify that Neumann BC hold
% u_y = diff(u,y,1);
% u_z = diff(u,z,1);
% 
% subs(u_y,y,'1') %==0
% subs(u_z,z,'1') %==0
           
%%
f = @(x,y,z) 2*x*exp(-(y - 1)^2*(z - 1)^2)*(y - 1)^2 + 2*x*exp(-(y - 1)^2*(z - 1)^2)*(z - 1)^2 - x*exp(-(y - 1)^2*(z - 1)^2)*(2*y - 2)^2*(z - 1)^4 - x*exp(-(y - 1)^2*(z - 1)^2)*(2*z - 2)^2*(y - 1)^4;

F = zeros(size(X));
for i = 1:size(X,1)
    x = X(i);
    y = Y(i);
    z = Z(i);
    F(i) = f(x,y,z);
end


end