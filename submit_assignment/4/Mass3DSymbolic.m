function M = Mass3DSymbolic()

xi = sym('xi', 'real');
eta = sym('eta', 'real');
zeta = sym('zeta', 'real');

dx = sym('dx', 'real');
dy = sym('dy', 'real');
dz = sym('dz', 'real');


c=[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
%jacobian of the transformation
J(1,1) = 0.5*dx;
J(2,2) = 0.5*dy;
J(3,3) = 0.5*dz;

for i=1:8
    N(i) = 1/8*( 1+c(i,1)*xi )*( 1+c(i,2)*eta )*( 1+c(i,3)*zeta );
end

F = det(J) * N' * N;
M = int(int(int(F, 'xi', -1, 1), 'eta', -1, 1), 'zeta', -1, 1);

end