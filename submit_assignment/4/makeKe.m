function Ke = makeKe(e, mesh)
% makeMe Constructs a local laplace matrix
% e identifies the element by its index in mesh.Elements array

dx = mesh.delta_x;
dy = mesh.delta_y;
dz = mesh.delta_z;

if (strcmp(mesh.grid_type,'hexahedra'))
    Ke = [ ... 
    (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz);
    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz);
    (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz);
    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz);
    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz);
    (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy);
  - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz);
    (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz)];

else
    display('Laplatian for tetrahedra elements not implemented');
%     %find coeff. matrix -> inverse of barycentric coordinates
%     nodes = mesh.Elements(e,:); 
%     coords = ones(3,3);
%     coords(:,1) = mesh.Points(nodes,1); %x coords
%     coords(:,2) = mesh.Points(nodes,2); %y coords
%     coeff = inv(coords);
%     
%     %make grad Ni; grad(Ni) * grad(Nj)
%     %make integral -> (ai*aj + bi*bj)*V
%     ai = coeff(1,:);
%     bi = coeff(2,:);
%     
%     %use outer product
%     Ke = ai'*ai + bi'*bi;
%     
%     %get volume of the element
%     Ve = 1/2 * abs(det(coords));
%     Ke = Ke * Ve;
end

end