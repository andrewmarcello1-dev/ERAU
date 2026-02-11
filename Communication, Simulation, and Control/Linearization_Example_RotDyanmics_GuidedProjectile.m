% Problem 1(d) phase-plane trajectories for linearized systems
% x1dot = 2*x1 - x1*x2
% x2dot = 2*x1^2 - x2

clc; clear; close all;

% Eq. Points

xe = [2;0;0];
ue = [0 0 0];

% State matrices

A = [0 0 0;
    0 0 -2;
    0 2 0];
dx01 = [0 1 0]';
dx02 = [0 0 -1]';

%% Linearized Function
exp_fun = @(t) [1 0 0;
             0 cos(2*t) -sin(2*t);
             0 sin(2*t) cos(2*t)];
f = @(t,dx0) exp_fun(t)*(dx0);


%% Exact Function
f_exact = @(t,x) [0;
                 -x(1)*x(3);
                 x(1)*x(2)];
u_fun = @(t,x) ue;

%% Calculations
tspan = [0 20];

x0_exact = [xe xe] + [dx01 dx02];

[t1, x_nl1] = ode45(f_exact, tspan, x0_exact(:,1));
[t2, x_nl2] = ode45(f_exact, tspan, x0_exact(:,2));

t = linspace(0,20,500)';
x_lin1 = zeros(length(t),3);
x_lin2 = x_lin1;

for k = 1:length(t)
    x_lin1(k,:) = xe + f(t(k), dx01);
    x_lin2(k,:) = xe + f(t(k), dx02);
end


% Plot results (actual state x(t))
figure; grid on; hold on;

subplot(2,1,1)
plot(t1, x_nl1(:,1), 'LineWidth', 2); hold on
plot(t1, x_nl1(:,2), 'LineWidth', 2); hold on
plot(t1, x_nl1(:,3), 'LineWidth', 2); hold on
plot(t, x_lin1(:,1), '--', 'LineWidth', 1); hold on
plot(t, x_lin1(:,2), '--', 'LineWidth', 1); hold on
plot(t, x_lin1(:,3), '--', 'LineWidth', 1);
xlabel('t (s)');
ylabel('\omega (rad/s)');
legend('\omega_x Exact','\omega_y Exact','\omega_z Exact','\omega_x Linearization','\omega_y Linearization','\omega_z Linearization');
title('Nonlinear (exact) rigid-body dynamics x_e = [0 1 0]^T');
xlim([0 20])
ylim([-2 3])

subplot(2,1,2)
plot(t2, x_nl2(:,1), 'LineWidth', 2); hold on
plot(t2, x_nl2(:,2), 'LineWidth', 2); hold on
plot(t2, x_nl2(:,3), 'LineWidth', 2); hold on
plot(t, x_lin2(:,1), '--', 'LineWidth', 1); hold on
plot(t, x_lin2(:,2), '--', 'LineWidth', 1); hold on
plot(t, x_lin2(:,3), '--', 'LineWidth', 1);
xlabel('t (s)');
ylabel('\omega (rad/s)');
legend('\omega_x','\omega_y','\omega_z');
title('Nonlinear (exact) rigid-body dynamics x_e = [0 0 -1]^T');
xlim([0 20])
ylim([-2 3])