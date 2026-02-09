%% ShockTube_2D.m
% This MATLAB script solves the 2D Euler equations for a shock tube 
% using a third‑order Runge–Kutta scheme with a finite-difference formulation
% based on Van Leer flux splitting. The computational domain is defined as
% x ∈ [-1,1] and y ∈ [–0.1,0.1] (yielding equal grid spacing) with 101 points in x 
% and 11 points in y. Inviscid wall (reflective) boundary conditions are applied
% on every wall.
%
% The conserved variables are:
%    U(:,:,1) = density (ρ)
%    U(:,:,2) = x‑momentum (ρ·u)
%    U(:,:,3) = y‑momentum (ρ·v)
%    U(:,:,4) = total energy (E = p/(γ–1) + 0.5·ρ·(u²+v²))
%
% The primitive variables are:
%    V(:,:,1) = density (ρ)
%    V(:,:,2) = x‑velocity (u)
%    V(:,:,3) = y‑velocity (v)
%    V(:,:,4) = pressure (p)

clear; close all; clc

%% Parameters

% Numerical parameters
CFL = 1;
nt = 50;                % Number of time steps
Nx = 101;               % Number of x-points
Ny = 11;                % Number of y-points

% Domain in x is [-1,1]
x = linspace(-1, 1, Nx)';      
dx = (x(end)-x(1))/(Nx-1);   % x-grid spacing

% For equal grid spacing, choose y-domain so that Δy = Δx.
dy = dx;                   % y-grid spacing equal to dx.
% Center the narrow y-domain around 0.
y_center = 0;
L_y = dy*(Ny-1);
y = linspace(y_center - L_y/2, y_center + L_y/2, Ny)';

% Gas properties
G  = 1.8;           % Ratio of specific heats (γ)
R  = 287;           % Gas constant (J/(kg·K))
Cv = R/(G-1);       % Constant-volume specific heat

%% Initial Conditions

% Left Region (x < 0): state 1
rho1 = 1.2;       % Density (kg/m^3)
p1   = 1e5;       % Pressure (Pa)
u1   = 0;         % x-velocity (m/s)
v1   = 0;         % y-velocity (m/s)

% Right Region (x >= 0): state 2
rho2 = 2.4;       % Density (kg/m^3)
p2   = 2e5;       % Pressure (Pa)
u2   = 0;         % x-velocity (m/s)
v2   = 0;         % y-velocity (m/s)

% Initialize the primitive variable array V (Nx x Ny x 4)
% V(:,:,1)=ρ, V(:,:,2)=u, V(:,:,3)=v, V(:,:,4)=p.
V_old = zeros(Nx, Ny, 4);
for i = 1:Nx
    for j = 1:Ny
        if x(i) < 0
            V_old(i,j,1) = rho1;
            V_old(i,j,4) = p1;
        else
            V_old(i,j,1) = rho2;
            V_old(i,j,4) = p2;
        end
        V_old(i,j,2) = 0;  % initial x-velocity
        V_old(i,j,3) = 0;  % initial y-velocity
    end
end

% Compute conserved variables U from V:
% U(:,:,1)=ρ, U(:,:,2)=ρ·u, U(:,:,3)=ρ·v,
% U(:,:,4)= Cv·p/R + 0.5·ρ·(u²+v²)
U_old = zeros(Nx, Ny, 4);
U_old(:,:,1) = V_old(:,:,1);
U_old(:,:,2) = V_old(:,:,1).*V_old(:,:,2);
U_old(:,:,3) = V_old(:,:,1).*V_old(:,:,3);
U_old(:,:,4) = Cv*V_old(:,:,4)/R + 0.5*V_old(:,:,1).*(V_old(:,:,2).^2 + V_old(:,:,3).^2);

%% Storage for Time-Stepping
dt_arr = zeros(nt,1);
t_plot = zeros(nt,1);
t_elap = 0;

% Normalization references from the left state
rho_ref = rho1;
p_ref   = p1;

% Create meshgrid for the contour plots (meshgrid returns arrays of size Ny×Nx)
[X, Y] = meshgrid(x, y);

%% Analytical Solution (McCormick)

x_a = [0 0.1040 0.1040 0.4155 0.4155 0.7401 0.8416 1]; % Position in Shock tube
x_a = (x_a*2)-1; % Rescale from [0,1] to [-1,1]
r_a = [1.2 1.2 1.4421 1.4421 1.9645 1.9645 2.4 2.4]; % Density (kg/m^3)
p_a = [1 1 1.4018 1.4018 1.4018 1.4018 2 2]*1e5; % Pressure (bar)

% Normalization references (left-side)
rho_ref = rho1;
p_ref = p1; % in Pa

%% Create Animation Figure with Contour and 1D Centerline Plots
figure('Name','2D Shock Tube Simulation');

% Top left subplot: Contour plot for normalized density
subplot(2,2,1);
hCont_density = contourf(X, Y, squeeze(V_old(:,:,1))'/rho_ref, 20);
colorbar;
title('Normalized Density (Contour)');
xlabel('x (m)'); ylabel('y (m)');

% Top right subplot: Contour plot for normalized pressure
subplot(2,2,2);
hCont_pressure = contourf(X, Y, squeeze(V_old(:,:,4))'/p_ref, 20);
colorbar;
title('Normalized Pressure (Contour)');
xlabel('x (m)'); ylabel('y (m)');

% For the centerline, choose the middle index in y.
j_center = round(Ny/2);

% Bottom left subplot: 1D centerline for normalized density
subplot(2,2,3);
hLine_density = plot(x, squeeze(V_old(:,j_center,1))/rho_ref, 'b-', 'LineWidth', 2);
title('Normalized Density (Centerline)');
xlabel('x (m)'); ylabel('\rho/\rho_0');
grid on;

% Bottom right subplot: 1D centerline for normalized pressure
subplot(2,2,4);
hLine_pressure = plot(x, squeeze(V_old(:,j_center,4))/p_ref, 'r-', 'LineWidth', 2);
title('Normalized Pressure (Centerline)');
xlabel('x (m)'); ylabel('p/p_0');
grid on;

Tmat = zeros(Nx,Ny); % Temperature for finding speed of sound later

%% RK3 Time Integration Loop

for z = 1:nt
    % Reconstruct primitive variables from U_old.
    rho = U_old(:,:,1);
    u   = U_old(:,:,2)./rho;
    v   = U_old(:,:,3)./rho;
    E   = U_old(:,:,4);
    p   = (G-1)*(E - 0.5*rho.*(u.^2+v.^2));
    
    % Compute local speed of sound and dt from CFL condition.
    c = sqrt(G * p ./ rho);
    % Use dt = CFL*dx / max(|u|+|v|+c) over all grid points.
    dt = CFL * dx / max(max(abs(u) + abs(v) + c));
    dt_arr(z) = dt;
    t_elap = t_elap + dt;
    t_plot(z) = t_elap;
    
    U_n = U_old;  % Save a copy for RK stages
    
    % RK3 stages
    for k = 1:3

        for i = 1:Nx
            for j = 1:Ny
                Tmat(i,j) = V_old(i,j,4)/(R*V_old(i,j,1)); % Resolve Temperature at every point
            end
        end
        alp = 1/(5 - k);  % Stage weights: 1/4, 1/3, 1/2 for k=1,2,3
        
        % Compute the 2D residual using the Van Leer FD scheme.
        Res = computeRHS2D(U_old, dx, dy, G);
        
        % Update the conserved variables for this stage.
        U_new = U_old - alp * dt * Res;
        
        % Recover primitive variables V_new from U_new.
        rho_new = U_new(:,:,1);
        u_new = U_new(:,:,2)./rho_new;
        v_new = U_new(:,:,3)./rho_new;
        p_new = (G-1)*(U_new(:,:,4) - 0.5*rho_new.*(u_new.^2+v_new.^2));
        V_new = cat(3, rho_new, u_new, v_new, p_new);
        
        % Apply Boundary Conditions (Inviscid Wall)

        % Velocity (x)
        V_new(:,1,2) = V_new(:,2,2);
        V_new(:,end,2) = V_new(:,end-1,2);
        V_new(1,:,2) = V_new(2,:,2);
        V_new(end,:,2) = V_new(end-1,:,2);

        % Velocity (y)
        V_new(:,1,3) = 0;
        V_new(:,end,3) = 0;
        V_new(1,:,3) = 0;
        V_new(end,:,3) = 0;

        % Pressure
        V_new(:,1,4) = V_new(:,2,4);
        V_new(:,end,4) = V_new(:,end-1,4);
        V_new(1,:,4) = V_new(2,:,4);
        V_new(end,:,4) = V_new(end-1,:,4);

        % Temperature
        Tmat(:,1) = Tmat(:,2);
        Tmat(1,:) = Tmat(2,:);
        Tmat(:,end) = Tmat(:,end-1);
        Tmat(end,:) = Tmat(end-1,:);

        % Density
        V_new(:,1,1) = V_new(:,1,4)./(R*Tmat(:,1));
        V_new(:,end,1) = V_new(:,end,4)./(R*Tmat(:,end));
        V_new(1,:,1) = V_new(1,:,4)./(R*Tmat(1,:));
        V_new(end,:,1) = V_new(end,:,4)./(R*Tmat(end,:));
        
        % Recompute the conserved variables from updated V_new.
        U_new(:,:,1) = V_new(:,:,1);
        U_new(:,:,2) = V_new(:,:,1).*V_new(:,:,2);
        U_new(:,:,3) = V_new(:,:,1).*V_new(:,:,3);
        U_new(:,:,4) = Cv*V_new(:,:,4)/R + 0.5*V_new(:,:,1).*(V_new(:,:,2).^2 + V_new(:,:,3).^2);
        
        % Update for next RK stage.
        U_old = U_new;
        V_old = V_new;
    end
    
    % Update the contour subplots.
    subplot(2,2,1);
    % Update density contour.
    cla;
    hCont_density = contourf(X, Y, squeeze(V_old(:,:,1))'/rho_ref, 20);
    colorbar;
    title(sprintf('Normalized Density (Contour) at t = %.4f s', t_elap));
    xlabel('x (m)'); ylabel('y (m)');
    
    subplot(2,2,2);
    % Update pressure contour.
    cla;
    hCont_pressure = contourf(X, Y, squeeze(V_old(:,:,4))'/p_ref, 20);
    colorbar;
    title(sprintf('Normalized Pressure (Contour) at t = %.4f s', t_elap));
    xlabel('x (m)'); ylabel('y (m)');
    
    % Update the centerline plots.
    % Extract centerline data (using the middle index in y).
    centerDensity = squeeze(V_old(:,j_center,1))/rho_ref;
    centerPressure = squeeze(V_old(:,j_center,4))/p_ref;
    
    subplot(2,2,3);
    set(hLine_density, 'YData', centerDensity);
    title(sprintf('Normalized Density (Centerline) at t = %.4f s', t_elap));
    
    subplot(2,2,4);
    set(hLine_pressure, 'YData', centerPressure);
    title(sprintf('Normalized Pressure (Centerline) at t = %.4f s', t_elap));
    
    drawnow;
    % pause(0.05); % Uncomment to slow down the animation if desired
end

%% Function: Compute 2D Residual via Van Leer FD Scheme
function Res = computeRHS2D(U, dx, dy, G)
    % U: (Nx x Ny x 4) array of conserved variables.
    [Nx, Ny, ~] = size(U);
    % Create an extended array with ghost cells (size: (Nx+2) x (Ny+2) x 4)
    U_ext = zeros(Nx+2, Ny+2, 4);
    U_ext(2:Nx+1, 2:Ny+1, :) = U;
    
    % Apply inviscid wall boundary conditions on every wall:
    % Left boundary (i=1)
    for j = 2:Ny+1
        U_ext(1,j,1) = U_ext(2,j,1);
        U_ext(1,j,2) = -U_ext(2,j,2);  % Reflect u
        U_ext(1,j,3) = U_ext(2,j,3);
        U_ext(1,j,4) = U_ext(2,j,4);
    end
    % Right boundary (i=Nx+2)
    for j = 2:Ny+1
        U_ext(Nx+2,j,1) = U_ext(Nx+1,j,1);
        U_ext(Nx+2,j,2) = -U_ext(Nx+1,j,2);
        U_ext(Nx+2,j,3) = U_ext(Nx+1,j,3);
        U_ext(Nx+2,j,4) = U_ext(Nx+1,j,4);
    end
    % Bottom boundary (j=1)
    for i = 2:Nx+1
        U_ext(i,1,1) = U_ext(i,2,1);
        U_ext(i,1,2) = U_ext(i,2,2);
        U_ext(i,1,3) = -U_ext(i,2,3);  % Reflect v
        U_ext(i,1,4) = U_ext(i,2,4);
    end
    % Top boundary (j=Ny+2)
    for i = 2:Nx+1
        U_ext(i,Ny+2,1) = U_ext(i,Ny+1,1);
        U_ext(i,Ny+2,2) = U_ext(i,Ny+1,2);
        U_ext(i,Ny+2,3) = -U_ext(i,Ny+1,3);
        U_ext(i,Ny+2,4) = U_ext(i,Ny+1,4);
    end
    % Corners (set as adjacent interior)
    U_ext(1,1,:) = U_ext(2,2,:);
    U_ext(1,Ny+2,:) = U_ext(2,Ny+1,:);
    U_ext(Nx+2,1,:) = U_ext(Nx+1,2,:);
    U_ext(Nx+2,Ny+2,:) = U_ext(Nx+1,Ny+1,:);
    
    % Compute primitive variables in the extended domain.
    rho_ext = U_ext(:,:,1);
    u_ext = U_ext(:,:,2)./rho_ext;
    v_ext = U_ext(:,:,3)./rho_ext;
    E_ext = U_ext(:,:,4);
    p_ext = (G-1)*(E_ext - 0.5*rho_ext.*(u_ext.^2 + v_ext.^2));
    
    % Compute the Van Leer flux splitting in the x-direction.
    [F_plus, F_minus] = van_leer_flux_x(rho_ext, u_ext, v_ext, p_ext, E_ext, G);
    % Compute the Van Leer flux splitting in the y-direction.
    [G_plus, G_minus] = van_leer_flux_y(rho_ext, u_ext, v_ext, p_ext, E_ext, G);
    
    % Compute the residual for interior cells (i = 2:Nx+1, j = 2:Ny+1).
    Res = zeros(Nx, Ny, 4);
    for i = 2:Nx+1
        for j = 2:Ny+1
            dF = ( F_plus(i,j,:) - F_plus(i-1,j,:) + F_minus(i+1,j,:) - F_minus(i,j,:) ) / dx;
            dG = ( G_plus(i,j,:) - G_plus(i,j-1,:) + G_minus(i,j+1,:) - G_minus(i,j,:) ) / dy;
            Res(i-1,j-1,:) = dF + dG;
        end
    end
end

%% Function: Van Leer Flux Splitting in x-direction
function [F_plus, F_minus] = van_leer_flux_x(rho, u, v, p, E, G)
    % Compute the local speed of sound and Mach number in the x-direction.
    c = sqrt(G * p ./ rho);
    M = u ./ c;
    
    % Initialize splitting coefficients.
    a_plus = zeros(size(M));
    a_minus = zeros(size(M));
    b_plus = zeros(size(M));
    b_minus = zeros(size(M));
    
    % For |M| <= 1.
    mask = abs(M) <= 1;
    a_plus(mask) = 0.25*(M(mask)+1).^2;
    a_minus(mask) = -0.25*(M(mask)-1).^2;
    b_plus(mask) = 0.25*(M(mask)+1).^2.*(2 - M(mask));
    b_minus(mask) = 0.25*(M(mask)-1).^2.*(2 + M(mask));
    
    % For M > 1.
    mask = M > 1;
    a_plus(mask) = M(mask);
    a_minus(mask) = 0;
    b_plus(mask) = 0.5*(M(mask)+1);
    b_minus(mask) = 0;
    
    % For M < -1.
    mask = M < -1;
    a_plus(mask) = 0;
    a_minus(mask) = M(mask);
    b_plus(mask) = 0;
    b_minus(mask) = 0.5*(1 - M(mask));
    
    % Compute the split fluxes in the x-direction.
    % The flux vector F = [ρ·u, ρ·u²+p, ρ·u·v, (E+p)·u].
    F_plus = zeros([size(rho), 4]);
    F_minus = zeros([size(rho), 4]);
    
    F_plus(:,:,1) = rho .* c .* a_plus;
    F_plus(:,:,2) = rho .* u .* c .* a_plus + p .* b_plus;
    F_plus(:,:,3) = rho .* v .* c .* a_plus;
    F_plus(:,:,4) = (E + p) .* c .* a_plus + p .* u .* b_plus;
    
    F_minus(:,:,1) = rho .* c .* a_minus;
    F_minus(:,:,2) = rho .* u .* c .* a_minus + p .* b_minus;
    F_minus(:,:,3) = rho .* v .* c .* a_minus;
    F_minus(:,:,4) = (E + p) .* c .* a_minus + p .* u .* b_minus;
end

%% Function: Van Leer Flux Splitting in y-direction
function [G_plus, G_minus] = van_leer_flux_y(rho, u, v, p, E, G)
    % Compute the local speed of sound and Mach number in the y-direction.
    c = sqrt(G * p ./ rho);
    M = v ./ c;
    
    % Initialize splitting coefficients.
    a_plus = zeros(size(M));
    a_minus = zeros(size(M));
    b_plus = zeros(size(M));
    b_minus = zeros(size(M));
    
    % For |M| <= 1.
    mask = abs(M) <= 1;
    a_plus(mask) = 0.25*(M(mask)+1).^2;
    a_minus(mask) = -0.25*(M(mask)-1).^2;
    b_plus(mask) = 0.25*(M(mask)+1).^2.*(2 - M(mask));
    b_minus(mask) = 0.25*(M(mask)-1).^2.*(2 + M(mask));
    
    % For M > 1.
    mask = M > 1;
    a_plus(mask) = M(mask);
    a_minus(mask) = 0;
    b_plus(mask) = 0.5*(M(mask)+1);
    b_minus(mask) = 0;
    
    % For M < -1.
    mask = M < -1;
    a_plus(mask) = 0;
    a_minus(mask) = M(mask);
    b_plus(mask) = 0;
    b_minus(mask) = 0.5*(1 - M(mask));
    
    % Compute the split fluxes in the y-direction.
    % The flux vector G = [ρ·v, ρ·u·v, ρ·v²+p, (E+p)·v].
    G_plus = zeros([size(rho), 4]);
    G_minus = zeros([size(rho), 4]);
    
    G_plus(:,:,1) = rho .* c .* a_plus;
    G_plus(:,:,2) = rho .* u .* c .* a_plus;
    G_plus(:,:,3) = rho .* v .* c .* a_plus + p .* b_plus;
    G_plus(:,:,4) = (E + p) .* c .* a_plus + p .* v .* b_plus;
    
    G_minus(:,:,1) = rho .* c .* a_minus;
    G_minus(:,:,2) = rho .* u .* c .* a_minus;
    G_minus(:,:,3) = rho .* v .* c .* a_minus + p .* b_minus;
    G_minus(:,:,4) = (E + p) .* c .* a_minus + p .* v .* b_minus;
end
