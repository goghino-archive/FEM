function Me = makeMe(e, mesh)
% makeMe Constructs a local element mass matrix
% e identifies the element by its index in mesh.Elements array

dx = mesh.delta_x;
dy = mesh.delta_y;

if (strcmp(mesh.grid_type,'quadrilaterals'))

    Me = ... 
    [  (dx*dy)/9, (dx*dy)/18, (dx*dy)/36, (dx*dy)/18;
       (dx*dy)/18,  (dx*dy)/9, (dx*dy)/18, (dx*dy)/36;
       (dx*dy)/36, (dx*dy)/18,  (dx*dy)/9, (dx*dy)/18;
       (dx*dy)/18, (dx*dy)/36, (dx*dy)/18,  (dx*dy)/9  ];

else
    %use formula for mass matrix
    Me = ...
        [1/6 1/12 1/12;
         1/12 1/6 1/12;
         1/12 1/12 1/6];
     
     %get volume of the element
     nodes = mesh.Elements(e,:);
     coords = ones(3,3);
     coords(:,1) = mesh.Points(nodes,1); %x coords
     coords(:,2) = mesh.Points(nodes,2); %y coords
     Ve = 1/2 * abs(det(coords));
     
     Me = Me * Ve;
end

end