function plotResult(mesh, u)
%should work only for quadrilateral grid

  nx = mesh.N_x;
  ny = mesh.N_y;
  lx = mesh.L_x;
  ly = mesh.L_y;
  [X,Y]=meshgrid(0:lx/nx:lx, 0:ly/ny:ly);
  temp=0;
  for i=1:nx+1,
      for j=1:ny+1,
        temp=temp+1;    
        U(i,j)=u(temp);
      end
  end
  %contourf(X,Y,U);
  figure
  surf(X,Y,U);
  zlim auto;
  title('solution')
  view(2), colorbar
end