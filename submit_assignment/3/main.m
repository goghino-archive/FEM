discretizations = [10 20 40 80];
error = zeros(size(discretizations,1),3);
i = 1;

for N = discretizations
    %% generate the grid
    mesh = makeGrid(1,1,N,N,'triangles');
    %mesh = makeGrid(1,1,N,N,'quadrilaterals');

    %% create FEM operators 
    [M, K, b] = assembleDiscreteOperators(mesh);
    K_noBoundary = K;

    %% impose boundary conditions for Dirichlet boundaries
    %e.g. point 10 is on dirichlet boundary
    %K(10,:) = 0
    %K(10,10) = 1 %diagonal
    %b(10) = u0(x10, y10);

    markers = reshape(mesh.PointMarkers,[mesh.N,1]);
    boundaryPoints = find(markers);

    K(boundaryPoints,:) = 0;
    for p = boundaryPoints'
        K(p,p) = 1; %diagonal
    end

    X = mesh.Points(boundaryPoints,1);
    Y = mesh.Points(boundaryPoints,2);
    b(boundaryPoints) = u0(X,Y);

    %% solve
    u_h = K\b;

    %% create exact solution
    X = mesh.Points(:,1);
    Y = mesh.Points(:,2);
    u = u0(X,Y);

    %% visualize solution
    %writeMeshToVTKFile(mesh, u_h, 'solution')
   
    %% Error norms
    error(i, 1) = sqrt( (u - u_h)' * (u - u_h) ); % euclidian norm
    error(i, 2) = sqrt( (u - u_h)' * M * (u - u_h) ); % L2 norm
    error(i, 3) = sqrt((u - u_h)' * M * (u - u_h) + (u - u_h)' * K_noBoundary * (u - u_h)); %H1 norm
    i = i+1;

end

loglog(discretizations,error,'+-')
legend('Euclidean','L2','H1');
xlabel('Discretization NxN')
ylabel('Error Norm')
set(gca,'fontsize',12)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);