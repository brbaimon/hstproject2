%% =========================================================
% HST PROJECT 2
% Observability, Observer Design, and Nonlinear Simulations
% =========================================================

clc
close all
clear all

load('system_model.mat');   % brings in A, B, C, D, x_star, u_star, F, etc.

n = size(A, 1);   % number of states

%% =========================================================
% STEP 1: OBSERVABILITY
% =========================================================

Ob = obsv(A, C);

fprintf('========================================\n');
fprintf('       OBSERVABILITY ANALYSIS\n');
fprintf('========================================\n');

fprintf('\nObservability matrix Ob =\n');
disp(Ob)

fprintf('Rank(Ob) = %d (n = %d)\n', rank(Ob), n);

if rank(Ob) == n
    fprintf('Result: System is fully observable.\n');
else
    fprintf('Result: System is NOT fully observable.\n');
    fprintf('Checking detectability instead (unobservable modes must be stable)...\n');
end

%% =========================================================
% STEP 2: OBSERVER DESIGN (Luenberger Observer)
% =========================================================

% Reference: teammate's closed-loop eigenvalues were
% [0.6420, 0.8614, 0.9018, 0.8994] (slowest = 0.9018).
% Observer poles chosen ~2-3x faster (smaller magnitude)
% than the slowest closed-loop pole, so state estimation
% settles well before the controller needs it.

desired_observer_poles = [0.2, 0.25, 0.3, 0.35];

% Observer gain via pole placement on the dual system (A', C')
L = place(A', C', desired_observer_poles)';

% Resulting observer error-dynamics matrix
A_observer = A - L*C;

observer_eig = eig(A_observer);

fprintf('\n========================================\n');
fprintf('       OBSERVER DESIGN (Luenberger)\n');
fprintf('========================================\n');

fprintf('\nReference — closed-loop control poles (from teammate):\n');
disp([0.6420, 0.8614, 0.9018, 0.8994])

fprintf('Desired observer poles (faster than control):\n');
disp(desired_observer_poles)

fprintf('Observer gain L =\n');
disp(L)

fprintf('A - L*C =\n');
disp(A_observer)

fprintf('Observer eigenvalues:\n');
disp(observer_eig)

if all(abs(observer_eig) < 1)
    fprintf('Result: Observer is asymptotically stable — state estimate converges.\n');
    fprintf('Convergence is faster than the closed-loop control dynamics.\n');
else
    fprintf('Result: Observer is NOT stable — check pole placement.\n');
end


%% =========================================================
% STEP 3: OBSERVER-BASED CLOSED-LOOP SIMULATION (LINEAR)
% =========================================================

% Load K from the state-feedback design (teammate's script
% should also save K to system_model.mat — see note below)

N = 50;   % number of time steps to simulate

% Initial conditions
x0     = x_star + [0.5; -0.3; 0.2; -0.1];   % true state, perturbed from equilibrium
xhat0  = x_star;                             % observer starts AT equilibrium (no info yet)

x    = zeros(n, N);
xhat = zeros(n, N);
y    = zeros(size(C,1), N);
u    = zeros(1, N);

x(:,1)    = x0;
xhat(:,1) = xhat0;

for t = 1:N-1
    % Measured output (deviation form, since A/B/C are linearised
    % around x_star)
    y(:,t) = C * (x(:,t) - x_star);

    % Control law uses the ESTIMATED state, not the true state
    u(t) = -K * (xhat(:,t) - x_star);

    % True system evolves (linear model, deviation form)
    x(:,t+1) = x_star + A*(x(:,t) - x_star) + B*u(t);

    % Observer update
    xhat(:,t+1) = x_star + A*(xhat(:,t) - x_star) + B*u(t) ...
        + L*(y(:,t) - C*(xhat(:,t) - x_star));
end

estimation_error = x - xhat;

%% Plot results
figure;
subplot(2,1,1)
plot(1:N, x(1,:), 'b-', 1:N, xhat(1,:), 'r--')
legend('True p1', 'Estimated p1')
title('True vs Estimated State: p1')
xlabel('Time step')

subplot(2,1,2)
plot(1:N, vecnorm(estimation_error))
title('Norm of Estimation Error ||x - xhat||')
xlabel('Time step')