%% =========================================================
% FULL PIPELINE: Observer + Controller + Nonlinear System
% (unstable regime, omega = 0.063)
% =========================================================

clc
close all
clear all

load('system_model.mat');   % A, B, C, D, x_star, F, K, A_sym, B_sym, etc.

N = 200;
x0 = x_star + [0.05; -0.05; 0.02; -0.02];

om_unstable = 0.063;

%% =========================================================
% RE-LINEARISE AT THE UNSTABLE REGIME AND DESIGN K, L FOR IT
% =========================================================

syms p1 p2 b1 b2 u omega1 omega2 mu1 mu2 real

A_unstable = double(subs(A_sym, ...
    [p1 p2 b1 b2 u omega1 omega2 mu1 mu2], ...
    [p1_star p2_star b1_star b2_star 0 om_unstable om_unstable om_unstable om_unstable]));

B_unstable = double(subs(B_sym, ...
    [p1 p2 b1 b2 u omega1 omega2 mu1 mu2], ...
    [p1_star p2_star b1_star b2_star 0 om_unstable om_unstable om_unstable om_unstable]));

Q = eye(4);
R = 100;
K_unstable = dlqr(A_unstable, B_unstable, Q, R);

desired_observer_poles = [0.2, 0.25, 0.3, 0.35];  % same design choice as before
L_unstable = place(A_unstable', C', desired_observer_poles)';

fprintf('\nObserver gain L_unstable =\n');
disp(L_unstable)
fprintf('Observer eigenvalues (A_unstable - L_unstable*C):\n');
disp(eig(A_unstable - L_unstable*C))

%% =========================================================
% NONLINEAR SIMULATION: true state via F, estimate via LINEAR
% observer, control via ESTIMATED state
% =========================================================

x_true = zeros(4, N);
x_hat  = zeros(4, N);
y_meas = zeros(2, N);
u_hist = zeros(1, N);

x_true(:,1) = x0;
x_hat(:,1)  = x_star;   % observer starts with no information

for t = 1:N-1
    % Measurement from TRUE nonlinear state
    y_meas(:,t) = C * (x_true(:,t) - x_star);

    % Controller uses ESTIMATED state only (realistic — we don't
    % have access to the true state in practice)
    u_hist(t) = -K_unstable * (x_hat(:,t) - x_star);

    % TRUE system evolves nonlinearly
    xt = x_true(:,t);
    x_true(:,t+1) = F(xt(1), xt(2), xt(3), xt(4), u_hist(t), ...
        om_unstable, om_unstable, om_unstable, om_unstable);

    % Observer update (LINEAR model — standard approximation)
    x_hat(:,t+1) = x_star + A_unstable*(x_hat(:,t) - x_star) ...
        + B_unstable*u_hist(t) ...
        + L_unstable*(y_meas(:,t) - C*(x_hat(:,t) - x_star));
end

estimation_error_ctrl = x_true - x_hat;

%% =========================================================
% PLOTS: full pipeline result — ALL STATES
% =========================================================

state_names = {'p1', 'p2', 'b1', 'b2'};

figure;
for i = 1:4
    subplot(2,2,i)
    plot(1:N, x_true(i,:), 'b-', 1:N, x_hat(i,:), 'r--')
    legend(['True ' state_names{i}], ['Estimated ' state_names{i}])
    title(['True vs Estimated: ' state_names{i}])
    xlabel('t')
end
sgtitle('Full Pipeline: True vs Estimated States under Observer-Based Control')
% saveas(gcf, 'fig5_pipeline_allstates.png')

figure;
plot(1:N, vecnorm(estimation_error_ctrl))
title('Estimation Error Norm ||x - xhat|| (all states)')
xlabel('t')
% saveas(gcf, 'fig6_pipeline_error.png')
