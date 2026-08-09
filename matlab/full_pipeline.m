%% =========================================================
% FULL PIPELINE: Observer + Controller + Nonlinear System
% (unstable regime, ω = 0.063)
% =========================================================

% Reuse A_unstable, B_unstable, K_unstable from before.
% Design an observer gain L_unstable for the SAME regime.

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

    % Observer update (LINEAR model — this is the standard
    % approximation; the observer doesn't know the true nonlinear f)
    x_hat(:,t+1) = x_star + A_unstable*(x_hat(:,t) - x_star) ...
        + B_unstable*u_hist(t) ...
        + L_unstable*(y_meas(:,t) - C*(x_hat(:,t) - x_star));
end

estimation_error_ctrl = x_true - x_hat;

%% =========================================================
% PLOTS: full pipeline result
% =========================================================

% figure;
% subplot(2,1,1)
% plot(1:N, x_true(1,:), 'b-', 1:N, x_hat(1,:), 'r--')
% legend('True p1 (nonlinear)', 'Estimated p1 (observer)')
% title('Full Pipeline: True vs Estimated p1 under Observer-Based Control')
% xlabel('t')
% 
% subplot(2,1,2)
% plot(1:N, vecnorm(estimation_error_ctrl))
% title('Estimation Error Norm ||x - xhat||')
% xlabel('t')

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

figure;
plot(1:N, vecnorm(estimation_error_ctrl))
title('Estimation Error Norm ||x - xhat|| (all states)')
xlabel('t')