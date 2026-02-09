% Problem 1(b) phase-plane trajectories for nonlinear system
% x1dot = 2*x1 - x1*x2
% x2dot = 2*x1^2 - x2

clc; clear; close all;

f = @(t,x) [ 2*x(1) - x(1)*x(2);      % x1dot
             2*x(1)^2 - x(2) ];      % x2dot

tspan = [0 10];

% Choose a grid of initial conditions (edit as desired)
x1_0 = linspace(-3,3,14);
x2_0 = linspace(-6,6,26);

figure; hold on; grid on;
xlabel('x_1'); ylabel('x_2');
title('Phase-plane trajectories: nonlinear system');

for i = 1:length(x1_0)
    for j = 1:length(x2_0)
        x0 = [x1_0(i); x2_0(j)];
        [t,x] = ode45(f, tspan, x0);
        plot(x(:,1), x(:,2), 'LineWidth', 1.0);
    end
end

% Mark equilibrium points
plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',6);
plot(1,2,'ko','MarkerFaceColor','k','MarkerSize',6);
plot(-1,2,'ko','MarkerFaceColor','k','MarkerSize',6);

% Optional: add direction arrows (quiver field)
[x1g,x2g] = meshgrid(linspace(-7,7,50), linspace(-6,10,50));
u = 2*x1g - x1g.*x2g;
v = 2*x1g.^2 - x2g;
quiver(x1g,x2g,u,v,0.7);

axis([-3 3 -6 6]);
