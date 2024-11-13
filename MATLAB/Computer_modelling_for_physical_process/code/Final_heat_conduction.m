clc;
clear all;
close all;
% Constants and parameters
Lx = 2;  %Length of the plate
Ly = 1;  % Width of the plate
thickness = 0.15;  
c = 50;  % Thermal conductivity 
Q1 = -200;  % Top side heat flux 
Q2 = 0;  % Left side heat flux 
T1 = 50;  % Right side temperature value
T2 = 10;  % Bottom side temperature value

Nx = 101;  %Number of grid points in the x-direction
Ny = 101;  %Number of grid points in the y-direction

% Grid spacing
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);

% Initializing the temperature matrix
T = 20 * ones(Ny, Nx); 

% Bottom side boundary condition
T(1, :) = T2;
    
% Top side boundary condition with heat flux
T(end, :) = T(end-1, :) - (Q1 * dy) / c;
    
% Right side boundary condition
T(:, end) = T1;
    
% Left side boundary condition with heat flux
T(:, 1) = T(:, 2) - (Q2 * dx) / c;

% Bottom right corner
T(end, end) = (T1 + T2)/2; 

% Applying boundary conditions inside the loop
max_iterations = 50000;
tolerance = 1e-5;

%Gauss-Seidel iterative method
for iter = 1:max_iterations
    T_old = T; 
    for i = 2:Ny
        for j = 1:Nx-1
          if(j==1)
                T(i,j)=T(i,j+1); %Left side boundary condition with heat flux 
          elseif(i==Ny)
              T(i,j)=T(i-1,j)-(Q1*dy)/c; %Top side boundary condition with heat flux
          else
            % Updating interior temperature values using the finite difference equation
            numer = ((dx)^2 * (T(i+1, j) + T_old(i-1, j))) + ((dy)^2 * (T_old(i, j+1) + T(i, j-1)));
            denom = 2 * (((dx)^2)+ ((dy)^2));
            T(i, j) = numer / denom;
            
          end 
        end
    end
    
 
    % Checking for convergence
    if max(abs(T - T_old), [], 'all') < tolerance
        break;
    end
end

% Plotting the temperature distribution with red-to-blue colormap
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
contourf(X, Y, T, 40, 'LineColor', 'none');
colormap(jet); 
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
title('Temperature Distribution in the Plate');
mid_horizontal = T(51,:); %to be changed based on the number of meshes
mid_verticle = T(:,51);   %to be changed based on the number of meshes