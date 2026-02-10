% Problem 1(d) phase-plane trajectories for linearized systems
% x1dot = 2*x1 - x1*x2
% x2dot = 2*x1^2 - x2

clc; clear; close all;

% Eq. Points

xe1 = [0;0];
xe2 = [1;2];
xe3 = [-1;2];

% State matrices

A1 = [2 0;
      0 -1];
A2 = [0 -1;
      4 -1];
A3 = [0 1;
     -4 -1];

f = @(t,x,A,xe) A*([x(1);x(2)]-xe);
f_exact = @(t,x) [ 2*x(1) - x(1)*x(2);
             2*x(1)^2 - x(2) ];

tspan = [0 10];

%Choose a grid of initial conditions (edit as desired)
x1_0 = linspace(-2,2,20);
x2_0 = linspace(-2,2,20);
x1_0_exact = linspace(-3,3,14);
x2_0_exact = linspace(-6,6,26);

%% Figure for eq1

figure; hold on; grid on;
xlabel('x_1'); ylabel('x_2');
title('Phase-plane trajectories: nonlinear system');

for i = 1:length(x1_0)
    for j = 1:length(x2_0)
        x0 = [x1_0(i); x2_0(j)];
        [t,x] = ode45(@(t,x) f(t,x,A1,xe1), tspan, x0);
        plot(x(:,1), x(:,2), 'LineWidth', 1.0);
    end
end

% Mark equilibrium points
plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',6);

xlim([-1e6 1e6])
ylim([-2.5e-3 2.5e-3])

%% Figure for eq2

figure; hold on; grid on;
xlabel('x_1'); ylabel('x_2');
title('Phase-plane trajectories: nonlinear system');

for i = 1:length(x1_0)
    for j = 1:length(x2_0)
        x0 = [x1_0(i); x2_0(j)];
        [t,x] = ode45(@(t,x) f(t,x,A2,xe2), tspan, x0);
        plot(x(:,1), x(:,2), 'LineWidth', 1.0);
    end
end

% Mark equilibrium points
plot(1,2,'ko','MarkerFaceColor','k','MarkerSize',6);

%% Figure for eq3

figure; hold on; grid on;
xlabel('x_1'); ylabel('x_2');
title('Phase-plane trajectories: nonlinear system');

for i = 1:length(x1_0)
    for j = 1:length(x2_0)
        x0 = [x1_0(i); x2_0(j)];
        [t,x] = ode45(@(t,x) f(t,x,A3,xe3), tspan, x0);
        plot(x(:,1), x(:,2), 'LineWidth', 1.0);
    end
end

% Mark equilibrium points
plot(-1,2,'ko','MarkerFaceColor','k','MarkerSize',6);

%% Figure for Exact

figure; hold on; grid on;
xlabel('x_1'); ylabel('x_2');
title('Phase-plane trajectories: nonlinear system');

for i = 1:length(x1_0_exact)
    for j = 1:length(x2_0_exact)
        x0 = [x1_0_exact(i); x2_0_exact(j)];
        [t,x] = ode45(f_exact, tspan, x0);
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