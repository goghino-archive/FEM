function be = makebe(e, mesh)
% makebe Constructs a function vector at the element e
% e identifies the element by its index in mesh.Elements array

I = mesh.Elements(e,:);
Xe = mesh.Points(I,1);
Ye = mesh.Points(I,2);
Ze = mesh.Points(I,3);

be = f(Xe,Ye,Ze);

end