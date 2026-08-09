%% =========================================================
% NONLINEAR SIMULATIONS
% Stable / Unstable / Chaotic / Controlled behaviour
% =========================================================

clc
close all
clear all

load('system_model.mat');   % A, B, C, D, x_star, F, K, A_sym, B_sym, etc.

N = 200;   % number of time steps

% Small perturbation from equilibrium as initial condition
x0 = x_star + [0.05; -0.05; 0.02; -0.02];

%% =========================================================
% REGIME DEFINITIONS (from bifurcation scan)
% =========================================================
%% =========================================================
% BIFURCATION SCAN (find stable / unstable / chaotic regimes)
% =========================================================

% clc
% close all
% clear all
% 
% load('system_model.mat');
% 
% N_transient = 300;   % steps to discard (let transient die out)
% N_keep      = 100;   % steps to keep and plot per omega value
% 
% omega_range = 0.01:0.005:0.20;   % sweep — adjust range if needed
% 
% bifurcation_p1 = [];   % will hold [omega, p1_value] pairs
% 
% for om = omega_range
% 
%     x = x_star + [0.05; -0.05; 0.02; -0.02];
%     diverged = false;
% 
%     for t = 1:N_transient
%         x = F(x(1), x(2), x(3), x(4), 0, om, om, om, om);
%         if any(abs(x) > 1e6) || any(isnan(x))
%             diverged = true;
%             break
%         end
%     end
% 
%     if diverged
%         continue   % skip this omega, too unstable to plot meaningfully
%     end
% 
%     for t = 1:N_keep
%         x = F(x(1), x(2), x(3), x(4), 0, om, om, om, om);
%         if any(abs(x) > 1e6) || any(isnan(x))
%             break
%         end
%         bifurcation_p1 = [bifurcation_p1; om, x(1)];
%     end
% end
% 
% figure;
% plot(bifurcation_p1(:,1), bifurcation_p1(:,2), '.', 'MarkerSize', 2)
% xlabel('\omega (speed parameter)')
% ylabel('p1 (long-run values)')
% title('Bifurcation diagram: p1 vs adjustment speed')

regimes = struct( ...
    'name',   {'Stable',  'Unstable', 'Chaotic'}, ...
    'omega1', {0.01,       0.063,      0.073}, ...
    'omega2', {0.01,       0.063,      0.073}, ...
    'mu1',    {0.01,       0.063,      0.073}, ...
    'mu2',    {0.01,       0.063,      0.073} ...
);

%% =========================================================
% RUN UNCONTROLLED NONLINEAR SIMULATIONS (3 regimes)
% =========================================================

results = struct();

for r = 1:numel(regimes)

    x = zeros(4, N);
    x(:,1) = x0;

    for t = 1:N-1
        xt = x(:,t);
        x(:,t+1) = F(xt(1), xt(2), xt(3), xt(4), 0, ...
                     regimes(r).omega1, regimes(r).omega2, ...
                     regimes(r).mu1,    regimes(r).mu2);
    end

    results.(regimes(r).name) = x;
end

%% =========================================================
% CONTROLLED CASE — controller designed FOR the unstable regime
% =========================================================

syms p1 p2 b1 b2 u omega1 omega2 mu1 mu2 real

om_unstable = 0.063;

A_unstable = double(subs(A_sym, ...
    [p1 p2 b1 b2 u omega1 omega2 mu1 mu2], ...
    [p1_star p2_star b1_star b2_star 0 om_unstable om_unstable om_unstable om_unstable]));

B_unstable = double(subs(B_sym, ...
    [p1 p2 b1 b2 u omega1 omega2 mu1 mu2], ...
    [p1_star p2_star b1_star b2_star 0 om_unstable om_unstable om_unstable om_unstable]));

Q = eye(4);
R = 100;
K_unstable = dlqr(A_unstable, B_unstable, Q, R);

x_ctrl = zeros(4, N);
x_ctrl(:,1) = x0;

for t = 1:N-1
    xt = x_ctrl(:,t);
    u_t = -K_unstable * (xt - x_star);
    x_ctrl(:,t+1) = F(xt(1), xt(2), xt(3), xt(4), u_t, ...
        om_unstable, om_unstable, om_unstable, om_unstable);
end

%% =========================================================
% LINEAR PREDICTION (for comparison, baseline speeds only)
% =========================================================

x_lin = zeros(4, N);
x_lin(:,1) = x0;

for t = 1:N-1
    x_lin(:,t+1) = x_star + A*(x_lin(:,t) - x_star);
end

%% =========================================================
% PLOTS
% =========================================================

figure;
subplot(2,2,1)
plot(1:N, results.Stable(1,:))
title('Stable regime: p1(t)')
xlabel('t')

subplot(2,2,2)
plot(1:N, results.Unstable(1,:))
title('Unstable regime: p1(t)')
xlabel('t')

subplot(2,2,3)
plot(1:N, results.Chaotic(1,:))
title('Chaotic regime: p1(t)')
xlabel('t')

subplot(2,2,4)
plot(1:N, x_ctrl(1,:))
title('Controlled (unstable speeds + regime-specific K): p1(t)')
xlabel('t')

figure;
plot(1:N, results.Stable(1,:), 'b-', 1:N, x_lin(1,:), 'r--')
legend('Nonlinear (stable regime)', 'Linear prediction')
title('Nonlinear vs Linear near equilibrium')
xlabel('t')