function F = makebe_vector(mesh)
% makebe_vector F is matrix N_v*N_e, each column is RHS for single element

F = zeros(mesh.N_v,mesh.N_e);

for e = 1:mesh.N_e
    I = mesh.Elements(e,:);
    Xe = mesh.Points(I,1);
    Ye = mesh.Points(I,2);
    Ze = mesh.Points(I,3);

    be = f(Xe,Ye,Ze);

    F(:,e) = be;
end

end