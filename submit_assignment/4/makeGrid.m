function[mesh] = makeGrid(L_x, L_y, L_z, N_x, N_y, N_z, grid_type)
% makeGrid Create a 3D grid for FEM
% Function generates a Nx - Ny - Nz grid of elements (not nodes/points). 
% The parameters L_x, L_y, L_z represent the length and height of the domain at each direction.
% 'triangles' generate a triangular grid with tetrahedra elements - NOT IMPLEMENTED
% 'quadrilaterals' generate a quadrilateral grid of hexahedra elements
%
% The output mesh is a MATLAB structure that contains all the information
% that are needed for the mesh generation and mesh description, namely:
% delta_x, delta_y, delta_z, N_x, N_y, N_z, L_x, L_y, L_z,
% N_e, N_v, N, Elements, Points, PointMarkers.

%%
% % OLD CODE - IT USES DIFFERENT INTERPRETATION OF N_X WHICH IS NUMBER OF
% % POINTS IN THE FINAL GRID, WHILE THE NEW VERSION BELOW INTERPRETS THIS
% % NUMBER AS NUBER OF ELEMENTS IN EACH DIM, THEREFOR NUMBER OF POINTS IS
% % N_X+1
%
% %length of segment between two nodes
% dx = L_x/(N_x-1);
% dy = L_y/(N_y-1);
% dz = L_z/(N_z-1);
% 
% a = 1;
% b = N_x;
% c = N_x*N_y;
% d = 0;
% 
% %number of vertices for different elements
% if(strcmp(grid_type, 'triangles'))
%     %N_v = 3;
%     %N_e = (N_x-1)*(N_y-1)*2;
%     display('Not Implemented for tetrahedra elements')
%     return;
% else %hexahedra
%     N_v = 8;
%     N_e = (N_x-1)*(N_y-1)*(N_z-1);
% end
% 
% %i goes in the x (column index) direction
% %j goes in the y (row index) direction
% %k goes in the z (horizontal slice index) direction
% nodes = zeros(N_x,N_y,N_z); %intersections of the grid lines
% coords = zeros(N_x*N_y*N_z, 3);%coordinates of nodes
% elements = zeros(N_e, N_v); %interior elements
% markersDiri = zeros(N_x*N_y*N_z,1); %flag for boundary points
% markersNeumann = zeros(N_x*N_y*N_z,1); %flag for boundary points
% id = 1;
% id_elem = 1;
% for k = 1:N_z
%     for j = 1:N_y
%         for i = 1:N_x
%             e = a*i + b*(j-1) + c*(k-1) + d;
%             nodes(j,i,k) = e;
%             coords(id,:) = [(i-1)*dx, (j-1)*dy, (k-1)*dz];
% 
%             %mark if Point belongs to the Dirichlet boundary
%             if (i==1 || j==1 || k==1 || i==N_x)
%                 markersDiri(id) = 1; %on boundary
%             else
%                 markersDiri(id) = 0;
%             end
%             
%             %mark if Point belongs to the Neumann boundary
%             %we don't need this
%             if (j==N_y || k==N_z )
%                 markersNeumann(id) = 1; %on boundary
%             else
%                 markersNeumann(id) = 0;
%             end
%             
%             id = id + 1;
% 
%             if (i<N_x && j<N_y && k<N_z)
%             if (strcmp(grid_type,'hexahedra'))
%                 %follow convenction when enumerating the corners
%                 %of the hexahedra element (VTK_HEXAHEDRON)
%                 elements(id_elem,:) = [e, e+1, e+N_x+1, e+N_x,...
%                     e+N_x*N_y, e+1+N_x*N_y, e+N_x+1+N_x*N_y, e+N_x+N_x*N_y];
%                 id_elem = id_elem + 1;
%             else %grid_type == 'tetrahedra'
%                 %elements(id_elem,:) = [e, e+1, e+N_x+1];
%                 %elements(id_elem + 1,:) = [e, e+N_x+1, e+N_x];
%                 %id_elem = id_elem + 2;
%             end
%             end
%         end
%     end
% end
% 
% 
% %delta_x, delta_y, N_x, N_y, L_x, L_y, N_e, N_v, N, Elements, Points, PointMarkers
% mesh = struct('delta_x', 0, 'delta_y', 0, 'delta_z', 0, 'N_x', 0, 'N_y', 0, 'N_z', 0, ...
%               'L_x', 0, 'L_y', 0, 'L_z', 0, 'N', 0, 'N_v', 0, 'N_e', 0, ...
%               'Points',[],'Elements',[],'grid_type','', 'PointMarkersDiri',[], 'PointMarkersNeumann',[]);
% 
% %construct the mesh structure
% if(not(N_x == N_y && N_y == N_z))
%     display('N_* dimensions not matching!');
%     return
% end
% 
% mesh.delta_x = dx;
% mesh.delta_y = dy;
% mesh.delta_z = dz;
% mesh.N_x = N_x;
% mesh.N_y = N_y;
% mesh.N_z = N_z;
% mesh.L_x = L_x;
% mesh.L_y = L_y;
% mesh.L_z = L_z;
% mesh.N_e = N_e;
% mesh.N_v = N_v;
% mesh.N = N_x*N_y*N_z;
% mesh.Elements = elements;
% mesh.Points = coords;
% mesh.grid_type = grid_type;
% mesh.PointMarkersDiri = markersDiri;
% mesh.PointMarkersNeumann = markersNeumann;

%%
mesh = struct('delta_x', 0, 'delta_y', 0, 'delta_z', 0, 'N_x', 0, 'N_y', 0, 'N_z', 0, ...
              'L_x', 0, 'L_y', 0, 'L_z', 0, 'N', 0, 'N_v', 0, 'N_e', 0, ...
              'Points',[],'Elements',[],'grid_type','', 'PointMarkersDiri',[], 'PointMarkersNeumann',[]);

%construct the mesh structure
if(not(N_x == N_y && N_y == N_z))
    display('N_* dimensions not matching!');
    return
end

%create grid of points
dx            = L_x/N_x;
dy            = L_y/N_y;
dz            = L_z/N_z;
x             = 0:dx:L_x;
y             = 0:dy:L_y;
z             = 0:dz:L_z;
[X, Y, Z]     = meshgrid(x, y, z); 
mesh.Points    = [ X(:), Y(:), Z(:) ];
Np = (N_x+1)*(N_y+1)*(N_z+1);
mesh.PointMarkers = zeros(Np, 1);
%mark if Point belongs to the Dirichlet boundary ( i==1 || i==N_x || j==1 || k==1 )
I = find(mesh.Points(:,1) == 0 | mesh.Points(:,1) == 1 | mesh.Points(:,2) == 0 | mesh.Points(:,3) == 0);
mesh.PointMarkersDiri(I) = 1;

%number of vertices for different elements
if(strcmp(grid_type, 'triangles'))
    %N_v = 3;
    %N_e = (N_x-1)*(N_y-1)*2;
    display('Not Implemented for tetrahedra elements')
    return;
else %hexahedra
    N_v = 8;
    N_e = (N_x)*(N_y)*(N_z);
end
mesh.grid_type = grid_type;

mesh.delta_x = dx;
mesh.delta_y = dy;
mesh.delta_z = dz;
mesh.N_x = N_x+1;
mesh.N_y = N_y+1;
mesh.N_z = N_z+1;
mesh.L_x = L_x;
mesh.L_y = L_y;
mesh.L_z = L_z;
mesh.N_e = N_e; %total number of elements
mesh.N_v = N_v; %number of vertices for each element
mesh.N = (N_x+1)*(N_y+1)*(N_z+1); %total number of nodes

%create element grid
k=1:N_z;
j=1:N_y;
i=1:N_x;

[I, J, K] = meshgrid(i, j, k);
I = I(:); J = J(:); K = K(:);

v1 = (K-1)*(N_x+1)*(N_y+1) + (J-1)*(N_x+1) + I;
v2 = (K-1)*(N_x+1)*(N_y+1) + (J-1)*(N_x+1) + I+1;
v3 = (K-1)*(N_x+1)*(N_y+1) +  J   *(N_x+1) + I+1;
v4 = (K-1)*(N_x+1)*(N_y+1) +  J   *(N_x+1) + I;

mesh.Elements     = [v1 v2 v3 v4  [v1 v2 v3 v4]+(N_x+1)*(N_y+1) ];

end