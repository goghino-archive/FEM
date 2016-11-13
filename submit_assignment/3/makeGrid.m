function[mesh] = makeGrid(L_x, L_y, N_x, N_y, grid_type)
% makeGrid Create a grid for FEM
% Function generates a Nx by Ny grid. 
% The parameters L_x, L_y represent the length and height of the domain at each direction.
% 'triangles' generate a triangular grid
% 'quadrilaterals' generate a quadrilateral grid
%
% The output mesh is a MATLAB structure that contains all the information
% that are needed for the mesh generation and mesh description, namely:
% delta_x, delta_y, N_x, N_y, L_x, L_y, N_e, N_v, N, Elements, Points, PointMarkers.

%length of segment between two nodes
dx = L_x/(N_x-1);
dy = L_y/(N_y-1);

a = 1;
b = N_x;
c = 0;

%number of vertices for different elements
if(strcmp(grid_type, 'triangles'))
    N_v = 3;
    N_e = (N_x-1)*(N_y-1)*2;
else %quadrilaterals
    N_v = 4;
    N_e = (N_x-1)*(N_y-1);
end

%i goes in the x (column index) direction
%j goes in the y (row index) direction
nodes = zeros(N_x,N_y); %intersections of the grid lines
coords = zeros(N_x*N_y, 2);%coordinates of nodes
elements = zeros(N_e, N_v); %interior elements
markers = zeros(N_x,N_y); %flag for boundary points
id_coord = 1;
id_elem = 1;
for j = 1:N_y
    for i = 1:N_x
        e = a*i + b*(j-1) + c;
        nodes(j,i) = e;
        coords(id_coord,:) = [(i-1)*dx, (j-1)*dy];
        id_coord = id_coord + 1;
        
        %mark if Point belongs to the boundary
        if (i==1 || j==1 || i==N_x || j==N_y)
            markers(j,i) = 1; %on boundary
        else
            markers(j,i) = 0;
        end
        
        if (i<N_x && j<N_y)
        if (strcmp(grid_type,'quadrilaterals'))
            %follow convenction when enumerating the corners of the element
            elements(id_elem,:) = [e, e+1, e+N_x+1, e+N_x];
            id_elem = id_elem + 1;
        else %grid_type == 'triangles'
            elements(id_elem,:) = [e, e+1, e+N_x+1];
            elements(id_elem + 1,:) = [e, e+N_x+1, e+N_x];
            id_elem = id_elem + 2;
        end
        end
    end
end

%construct the mesh structure
%delta_x, delta_y, N_x, N_y, L_x, L_y, N_e, N_v, N, Elements, Points, PointMarkers
mesh = struct('delta_x', 0, 'delta_y', 0, 'N_x', 0, 'N_y', 0, 'L_x', 0, 'L_y', 0, ...
'N', 0, 'N_v', 0, 'N_e', 0, 'Points',[],'Elements',[],'grid_type','', 'PointMarkers',[]);

mesh.delta_x = dx;
mesh.delta_y = dy;
mesh.N_x = N_x;
mesh.N_y = N_y;
mesh.L_x = L_x;
mesh.L_y = L_y;
mesh.N_e = N_e;
mesh.N_v = N_v;
mesh.N = N_x*N_y;
mesh.Elements = elements;
mesh.Points = coords;
mesh.grid_type = grid_type;
mesh.PointMarkers = markers;

%writeMeshToVTKFile(mesh, 'vtkfile');
end