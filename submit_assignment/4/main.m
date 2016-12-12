%find the analytical solution of f given u in Poisson equation:
% - \delta u(x,y,z) = f(x,y,z)
% for Dirichlet and Neumann BC:
% \partial u / \partial n = 0 for y=1, z=1
% u(x,y,z) = u0(x,y,z) elsewhere
% GOTO -> f.m
clear all;

discretizations = [10 15];
error = zeros(size(discretizations,1),3);
times = zeros(size(discretizations,1),2);
iter = 1;

for N = discretizations
    %% generate 3D hexa-hedron grid
    mesh = makeGrid(1,1,1,N,N,N,'hexahedra');

    %% create FEM operators 
    tstart = tic;
    %[M, K, b] = assembleDiscreteOperators(mesh);
    [M,K,b] = assembleFast(mesh);
    times(iter,1) = toc(tstart);
    
    K_noBoundary = K;

    %% impose boundary conditions for Dirichlet boundaries
    %e.g. point 10 is on dirichlet boundary
    %K(10,:) = 0
    %K(10,10) = 1 %diagonal
    %b(10) = u0(x10, y10);

    I = find(mesh.PointMarkersDiri);

    for i = 1:length(I);
      idx = I(i);    

      X = mesh.Points(idx,1);
      Y = mesh.Points(idx,2);
      Z = mesh.Points(idx,3);
    
      %fix matrix to keep symetry after applying dirichlet BC
      %that is, move it to RHS
      b = b - K(:,idx)*u0(X,Y,Z);
      b(idx) = u0(X,Y,Z);
      
      K(idx,:) = 0;
      K(:,idx) = 0;
      K(idx,idx) = 1.0;
    end
    
    
    %% impose boundary conditions for Neumann boundaries
    % it means adding some term to RHS (implied by weak formulation, where
    % Neumann BC shows up) but since it is zero in this specific case
    % it does not change RHS

    %% solve using iterative method (instead of direct solve u_h = K\b)
    tstart = tic;
    u_h = pcg(K,b,1e-10,1000);
    times(iter,2) = toc(tstart);

    %% create exact solution
    X = mesh.Points(:,1);
    Y = mesh.Points(:,2);
    Z = mesh.Points(:,3);
    u = u0(X,Y,Z);

    %% visualize solution
    %writeMeshToVTKFile(mesh, u, 'solution');
    %writeMeshToVTKFile(mesh, u_h, 'FEMsolution');
   
    %% Error norms
    error(iter, 1) = sqrt( (u - u_h)' * (u - u_h) ); % euclidian norm
    error(iter, 2) = sqrt( (u - u_h)' * M * (u - u_h) ); % L2 norm
    fprintf('L2 norm %f for grid size %dx%dx%d\n', error(iter, 2), N, N, N);
    error(iter, 3) = sqrt((u - u_h)' * M * (u - u_h) + (u - u_h)' * K_noBoundary * (u - u_h)); %H1 norm
    fprintf('H1 norm %f for grid size %dx%dx%d\n\n', error(iter, 3), N, N, N);
    iter = iter+1;

end

loglog(discretizations,error,'+-')
legend('Euclidean','L2','H1');
xlabel('Discretization NxN')
ylabel('Error Norm')

figure();
loglog(discretizations,times,'+-')
legend('Assembly','Solution');
xlabel('Discretization NxN')
ylabel('Time [s]')

set(gca,'fontsize',12)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
