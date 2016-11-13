function I = Laplace2DSymbolic()
xi = sym('xi', 'real');
eta = sym('eta', 'real');
dx = sym('dx', 'real');
dy = sym('dy', 'real');
c=[-1 -1; 1 -1; 1 1; -1 1];

J(1,1) = 0.5*dx;
J(2,2) = 0.5*dy;

for i=1:4
    N(i) = 1/4*( 1+c(i,1)*xi )*( 1+c(i,2)*eta );
end
Nx = diff(N, 'xi');
Ny = diff(N, 'eta');
dN = [Nx; Ny];
Jd = inv(J) * dN;
F = det(J)*Jd'*Jd;
I = int(int(F, 'xi', -1, 1), 'eta', -1, 1);

end