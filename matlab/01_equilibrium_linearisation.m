%% HST Project 2
% System Model, Equilibrium and Linearisation

clear;
clc;
close all;

%% 1. Parameters

beta = 1;
g1 = 1;
g2 = 1;
Am = 1;

eta1 = 0.1;
eta2 = 0.1;

r = 0.1;
m = 0.1;

epsilon1 = 0.5;

c1 = 0.1;
c2 = 0.1;

omega1 = 0.01;
omega2 = 0.01;

mu1 = 0.01;
mu2 = 0.01;

k0 = 0;

% Market-size parameter still needs to be verified
a = NaN;

%% 2. State variables
% x1 = p1
% x2 = p2
% x3 = b1
% x4 = b2

%% 3. Demand equations

% q1 = a - (1+r+m)*p1 + beta*(1+r)*p2 + g1*b1;
% q2 = a - (1+r)*p2 + beta*(1+r+m)*p1 + g2*b2;

%% 4. Profit equations

% pi1 = ((1+r+m)*p1 - c1)*q1 ...
%       - (eta1/2)*b1^2 - epsilon1*Am;

% pi2 = ((1+r)*p2 - c2)*q2 ...
%       - (eta2/2)*b2^2 - (1-epsilon1)*Am;

%% 5. Marginal-profit derivatives

%% 6. Nonlinear state equations

%% 7. Interior equilibrium

%% 8. Equilibrium verification

%% 9. Linearisation

%% 10. Display A, B and C matrices
