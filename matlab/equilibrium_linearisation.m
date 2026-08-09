%% =========================================================
% HST PROJECT 2
% System Model, Equilibrium and Linearisation
%
% State vector:
%
%       x = [p1; p2; b1; b2]
%
% p1 = Channel 1 retail price
% p2 = Channel 2 retail price
% b1 = Channel 1 blockchain maintenance effort
% b2 = Channel 2 blockchain maintenance effort
%
% Control input:
% u = additive price correction applied directly to p2
%
% Measured outputs:
% y = [p1; p2]
%% =========================================================

clear;
clc;
close all;


%% =========================================================
% STEP 1: PARAMETER VALUES
% ==========================================================

a = 10;                 % Initial market size

beta = 1;               % Competition intensity

g1 = 1;                 % Blockchain demand sensitivity, SC1
g2 = 1;                 % Blockchain demand sensitivity, SC2

Am = 1;                 % Fixed blockchain construction cost

% Calibrated blockchain maintenance-cost coefficients
% Paper states eta = 0.1, but this gives an economically
% infeasible negative interior equilibrium.
% eta = 100 gives a positive equilibrium consistent with
% the intended numerical behaviour.
eta1 = 100;
eta2 = 100;

r = 0.1;                % Tariff rate
m = 0.1;                % Commission rate

epsilon1 = 0.5;         % Blockchain fixed-cost sharing coefficient

c1 = 0.1;               % Unit cost, SC1
c2 = 0.1;               % Unit cost, SC2

omega1 = 0.01;          % Price adjustment speed, SC1
omega2 = 0.01;          % Price adjustment speed, SC2

mu1 = 0.01;             % Blockchain effort adjustment speed, SC1
mu2 = 0.01;             % Blockchain effort adjustment speed, SC2

% Fixed nominal delayed-feedback factor
% k0 = 0 gives the original baseline adjustment speeds.
k0 = 0;


%% =========================================================
% STEP 2: DEFINE SYMBOLIC VARIABLES
% ==========================================================

syms p1 p2 b1 b2 u real

% State vector
x = [p1;
     p2;
     b1;
     b2];


%% =========================================================
% STEP 3: DEMAND EQUATIONS
% ==========================================================

% Demand for Channel 1

q1 = a ...
     - (1+r+m)*p1 ...
     + beta*(1+r)*p2 ...
     + g1*b1;


% Demand for Channel 2

q2 = a ...
     - (1+r)*p2 ...
     + beta*(1+r+m)*p1 ...
     + g2*b2;


disp('========================================')
disp('DEMAND EQUATIONS')
disp('========================================')

disp('q1 = ')
pretty(q1)

disp('q2 = ')
pretty(q2)


%% =========================================================
% STEP 4: PROFIT EQUATIONS
% ==========================================================

% Channel 1 profit

pi1 = ((1+r+m)*p1 - c1)*q1 ...
      - (eta1/2)*b1^2 ...
      - epsilon1*Am;


% Channel 2 profit

pi2 = ((1+r)*p2 - c2)*q2 ...
      - (eta2/2)*b2^2 ...
      - (1-epsilon1)*Am;


disp('========================================')
disp('PROFIT EQUATIONS')
disp('========================================')

disp('pi1 = ')
pretty(expand(pi1))

disp('pi2 = ')
pretty(expand(pi2))


%% =========================================================
% STEP 5: MARGINAL PROFIT DERIVATIVES
% ==========================================================

% Price marginal profits

dpi1_dp1 = simplify(diff(pi1,p1));

dpi2_dp2 = simplify(diff(pi2,p2));


% Blockchain-effort marginal profits

dpi1_db1 = simplify(diff(pi1,b1));

dpi2_db2 = simplify(diff(pi2,b2));


disp('========================================')
disp('MARGINAL PROFIT DERIVATIVES')
disp('========================================')

disp('d(pi1)/d(p1) = ')
pretty(dpi1_dp1)

disp('d(pi2)/d(p2) = ')
pretty(dpi2_dp2)

disp('d(pi1)/d(b1) = ')
pretty(dpi1_db1)

disp('d(pi2)/d(b2) = ')
pretty(dpi2_db2)


%% =========================================================
% STEP 6: NONLINEAR STATE UPDATE EQUATIONS
% ==========================================================
%
% General form:
%
%       x(t+1) = f(x(t),u(t))
%
% The fixed delayed-feedback factor k0 scales the
% adjustment speeds.
%
% The group-defined control input u is added directly
% to the p2 update equation.
%% =========================================================


% p1 update

p1_next = p1 ...
          + (omega1/(1+k0))*p1*dpi1_dp1;


% p2 update + control input

p2_next = p2 ...
          + (omega2/(1+k0))*p2*dpi2_dp2 ...
          + u;


% b1 update

b1_next = b1 ...
          + (mu1/(1+k0))*b1*dpi1_db1;


% b2 update

b2_next = b2 ...
          + (mu2/(1+k0))*b2*dpi2_db2;


% Complete nonlinear state-space function

f = [p1_next;
     p2_next;
     b1_next;
     b2_next];


disp('========================================')
disp('NONLINEAR STATE UPDATE f(x,u)')
disp('========================================')

disp('p1(t+1) = ')
pretty(p1_next)

disp('p2(t+1) = ')
pretty(p2_next)

disp('b1(t+1) = ')
pretty(b1_next)

disp('b2(t+1) = ')
pretty(b2_next)


%% =========================================================
% STEP 7: SOLVE FOR POSITIVE INTERIOR EQUILIBRIUM
% ==========================================================
%
% At an interior equilibrium:
%
%       d(pi1)/d(p1) = 0
%       d(pi2)/d(p2) = 0
%       d(pi1)/d(b1) = 0
%       d(pi2)/d(b2) = 0
%
% and u* = 0.
%% =========================================================


sol = solve( ...
    dpi1_dp1 == 0, ...
    dpi2_dp2 == 0, ...
    dpi1_db1 == 0, ...
    dpi2_db2 == 0, ...
    [p1 p2 b1 b2]);


% Convert symbolic solution to numerical values

p1_star = double(sol.p1);

p2_star = double(sol.p2);

b1_star = double(sol.b1);

b2_star = double(sol.b2);


% Equilibrium state vector

x_star = [p1_star;
          p2_star;
          b1_star;
          b2_star];


%% =========================================================
% STEP 8: EQUILIBRIUM DEMANDS
% ==========================================================

q1_star = double(subs(q1, ...
    [p1 p2 b1 b2], ...
    [p1_star p2_star b1_star b2_star]));


q2_star = double(subs(q2, ...
    [p1 p2 b1 b2], ...
    [p1_star p2_star b1_star b2_star]));


%% =========================================================
% STEP 9: VERIFY EQUILIBRIUM
% ==========================================================

% Equilibrium control input

u_star = 0;


% Evaluate nonlinear system at x*

f_star = double(subs(f, ...
    [p1 p2 b1 b2 u], ...
    [p1_star p2_star b1_star b2_star u_star]));


% Fixed-point residual

equilibrium_residual = norm(f_star - x_star,2);


% Economic feasibility check

all_positive = all([ ...
    p1_star ...
    p2_star ...
    b1_star ...
    b2_star ...
    q1_star ...
    q2_star] > 0);


%% =========================================================
% DISPLAY EQUILIBRIUM RESULTS
% ==========================================================

fprintf('\n');

fprintf('========================================\n');
fprintf('       INTERIOR EQUILIBRIUM x*\n');
fprintf('========================================\n');

fprintf('p1* = %.10f\n',p1_star);
fprintf('p2* = %.10f\n',p2_star);
fprintf('b1* = %.10f\n',b1_star);
fprintf('b2* = %.10f\n',b2_star);


fprintf('\n');

fprintf('========================================\n');
fprintf('       DEMAND AT EQUILIBRIUM\n');
fprintf('========================================\n');

fprintf('q1* = %.10f\n',q1_star);
fprintf('q2* = %.10f\n',q2_star);


fprintf('\n');

fprintf('========================================\n');
fprintf('       FIXED-POINT VERIFICATION\n');
fprintf('========================================\n');

fprintf('||f(x*) - x*|| = %.12e\n', ...
    equilibrium_residual);


if equilibrium_residual < 1e-8

    fprintf('PASS: Equilibrium residual is below 1e-8.\n');

else

    fprintf('FAIL: Equilibrium residual is NOT below 1e-8.\n');

end


fprintf('\n');

fprintf('========================================\n');
fprintf('       ECONOMIC FEASIBILITY\n');
fprintf('========================================\n');

fprintf('p1* = %.6f\n',p1_star);
fprintf('p2* = %.6f\n',p2_star);
fprintf('b1* = %.6f\n',b1_star);
fprintf('b2* = %.6f\n',b2_star);
fprintf('q1* = %.6f\n',q1_star);
fprintf('q2* = %.6f\n',q2_star);


if all_positive

    fprintf('\nPASS: All states and demands are positive.\n');

else

    fprintf('\nWARNING: Equilibrium is NOT economically feasible.\n');

end


%% =========================================================
% STEP 10: LINEARISATION
% ==========================================================
%
% Around:
%
%       x = x*
%       u = u* = 0
%
% Define deviation variables:
%
%       Delta x = x - x*
%       Delta u = u - u*
%
% Linearised model:
%
% Delta x(t+1) = A Delta x(t) + B Delta u(t)
%
% where:
%
%       A = df/dx | x=x*
%       B = df/du | x=x*
%
%% =========================================================


%% 10.1 SYMBOLIC JACOBIAN A

A_sym = simplify(jacobian(f,x));


disp('========================================')
disp('SYMBOLIC JACOBIAN A(x)')
disp('========================================')

pretty(A_sym)


%% 10.2 NUMERICAL A MATRIX AT EQUILIBRIUM

A = double(subs(A_sym, ...
    [p1 p2 b1 b2 u], ...
    [p1_star p2_star b1_star b2_star 0]));


%% =========================================================
% STEP 11: INPUT MATRIX B
% ==========================================================
%
% u enters directly into the p2 update.
%
% Therefore:
%
%             [0]
%             [1]
%       B  =  [0]
%             [0]
%% =========================================================


B_sym = simplify(jacobian(f,u));


B = double(subs(B_sym, ...
    [p1 p2 b1 b2 u], ...
    [p1_star p2_star b1_star b2_star 0]));


%% =========================================================
% STEP 12: OUTPUT MATRIX C
% ==========================================================
%
% Only prices p1 and p2 are measured.
%
%       y = Cx
%
%% =========================================================


C = [1 0 0 0;
     0 1 0 0];


%% =========================================================
% STEP 13: DIRECT FEEDTHROUGH MATRIX D
% ==========================================================

% Output equation does not directly depend on u.

D = zeros(2,1);


%% =========================================================
% STEP 14: DISPLAY LINEARISED MATRICES
% ==========================================================

fprintf('\n');

fprintf('========================================\n');
fprintf('       LINEARISED STATE-SPACE MODEL\n');
fprintf('========================================\n');


fprintf('\nA MATRIX:\n');

disp(A)


fprintf('B MATRIX:\n');

disp(B)


fprintf('C MATRIX:\n');

disp(C)


fprintf('D MATRIX:\n');

disp(D)


fprintf('\n');

fprintf('Linearised system:\n\n');

fprintf('Delta x(t+1) = A Delta x(t) + B Delta u(t)\n');

fprintf('Delta y(t)   = C Delta x(t) + D Delta u(t)\n');


%% =========================================================
% STEP 15: FINAL SUMMARY
% ==========================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('              FINAL SUMMARY\n');
fprintf('========================================\n');

fprintf('\nEquilibrium state:\n');

disp(x_star)


fprintf('Equilibrium demands:\n');

fprintf('q1* = %.6f\n',q1_star);
fprintf('q2* = %.6f\n',q2_star);


fprintf('\nEquilibrium residual:\n');

fprintf('%.12e\n',equilibrium_residual);


fprintf('\nA =\n');

disp(A)


fprintf('B =\n');

disp(B)


fprintf('C =\n');

disp(C)


fprintf('========================================\n');
fprintf('        END OF YOUR ASSIGNED PART\n');
fprintf('========================================\n');

%% =========================================================
% STEP 16: SAVE RESULTS FOR DOWNSTREAM SCRIPTS
% ==========================================================

% Convert symbolic nonlinear update f(x,u) into a numeric function
% handle, so your nonlinear simulation script can call it directly
% without needing Symbolic Toolbox at runtime.
F = matlabFunction(f, 'Vars', [p1 p2 b1 b2 u]);

save('system_model.mat', ...
    'A', 'B', 'C', 'D', ...
    'x_star', 'u_star', ...
    'p1_star', 'p2_star', 'b1_star', 'b2_star', ...
    'q1_star', 'q2_star', ...
    'F');

fprintf('\nSaved system_model.mat with A, B, C, D, equilibrium values, and F (nonlinear map).\n');