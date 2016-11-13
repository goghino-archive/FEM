function[nodes,elements] = makeQuadMesh(nx,ny)

%number of nodes in each dimension
%nx = 100;
%ny = 100;

%length of segment between two nodes
dx = 1/(nx-1);
dy = 1/(ny-1);

a = 1;
b = nx;
c = 0;

%i goes in the x direction
%j goes in the y direction
nodes = zeros(nx,ny); %intersections of the grid lines
coords = zeros(nx*ny, 2);%coordinates of nodes
elements = zeros((nx-1)*(ny-1), 4); %interior elements
id_coord = 1;
id_elem = 1;
for j = 1:ny
    for i = 1:nx
        e = a*i + b*(j-1) + c;
        nodes(j,i) = e;
        coords(id_coord,:) = [(i-1)*dx, (j-1)*dy];
        id_coord = id_coord + 1;
        
        if (i<nx && j<ny)
            %follow convenction when enumerating the corners of the element
            elements(id_elem,:) = [e, e+1, e+nx+1, e+nx];
            id_elem = id_elem + 1;
        end
    end
end

%construct the mesh structure
myMesh = struct('Points',[],'Elements',[],'N',0,'N_v',0,'N_e',0,'grid_type','');
myMesh.N = nx*ny;
myMesh.Points = coords;
myMesh.N_v = 4;
myMesh.N_e = (nx-1)*(ny-1);
myMesh.Elements = elements;
myMesh.grid_type = 'quadrilaterals';

writeMeshToVTKFile(myMesh, 'vtkfile');

end