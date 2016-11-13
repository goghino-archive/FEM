function [Ms,Ks,b] = assembleFast(mesh)
dummy = 0;
Me = makeMe(dummy, mesh);
Ke = makeKe(dummy, mesh);
me = Me (:);
ke = Ke (:);
M = repmat(me, mesh.N_e, 1);
K = repmat(ke, mesh.N_e, 1);
i=1:mesh.N_v;
% ig = [1..N_v, 1..N_v, ... 1..N_v] N_v times
ig=repmat(i, 1, mesh.N_v);
% jg = [1 1 ... 1, 2 2 ... 2, ... , N_v N_v ... N_v]
jg=repmat(1:mesh.N_v, mesh.N_v, 1);
jg=jg(:)';
iA = mesh.Elements(:,ig)';
jA = mesh.Elements(:,jg)';
Ks = sparse(iA(:), jA(:), K, mesh.N, mesh.N);
Ms = sparse(iA(:), jA(:), M, mesh.N, mesh.N);

%F is matrix N_v*N_e, each column is RHS for single element
F = makebe_vector(mesh);
b = Me*F;

I = mesh.Elements';
I = I(:);
J = ones(mesh.N_v*mesh.N_e, 1);
B = sparse(I, J, b(:), mesh.N, 1);
b = full(B);

end