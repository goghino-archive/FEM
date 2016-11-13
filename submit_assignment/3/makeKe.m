function Ke = makeKe(e, mesh)
% makeMe Constructs a local laplace matrix
% e identifies the element by its index in mesh.Elements array

dx = mesh.delta_x;
dy = mesh.delta_y;

if (strcmp(mesh.grid_type,'quadrilaterals'))
    Ke = ... 
    [  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx), -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy);
        dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy),    dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy);
     -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy),  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx);
        dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy),    dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy)];

else
    %find coeff. matrix -> inverse of barycentric coordinates
    nodes = mesh.Elements(e,:); 
    coords = ones(3,3);
    coords(:,1) = mesh.Points(nodes,1); %x coords
    coords(:,2) = mesh.Points(nodes,2); %y coords
    coeff = inv(coords);
    
    %make grad Ni; grad(Ni) * grad(Nj)
    %make integral -> (ai*aj + bi*bj)*V
    ai = coeff(1,:);
    bi = coeff(2,:);
    
    %use outer product
    Ke = ai'*ai + bi'*bi;
    
    %get volume of the element
    Ve = 1/2 * abs(det(coords));
    Ke = Ke * Ve;
end

end