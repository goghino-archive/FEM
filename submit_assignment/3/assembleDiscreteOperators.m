function [M, K, b] = assembleDiscreteOperators(mesh)
% assembleDiscreteOperators
% Output
% M global mass matrix
% K global Laplacian
% b rhs multiplied by mass matrix

N = mesh.N;
N_e = mesh.N_e;
N_v = mesh.N_v;

M = zeros(N,N) ;
K = zeros(N,N) ;
b = zeros(N,1);

for e = 1:N_e
    % M_e element mass matrix
    % K_e element Laplacian
    % b_e element rhs
    Me = makeMe(e, mesh);
    Ke = makeKe(e, mesh);
    fp = makebe(e, mesh);
    for i = 1:N_v
        I = mesh.Elements(e, i);
        for j = 1:N_v
            J = mesh.Elements(e, j);
            M(I, J) = M(I, J) + Me(i, j);
            K(I, J) = K(I, J) + Ke(i, j);
        end
        be = Me*fp;
        b(I) = b(I) + be(i);
    end
end

end