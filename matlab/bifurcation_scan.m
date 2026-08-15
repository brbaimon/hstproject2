%% =========================================================
% BIFURCATION SCAN (find stable / unstable / chaotic regimes)
% =========================================================

clc
close all
clear all

load('system_model.mat');

N_transient = 300;   % steps to discard (let transient die out)
N_keep      = 100;   % steps to keep and plot per omega value

omega_range = 0.01:0.005:0.20;   % sweep — adjust range if needed

bifurcation_p1 = [];   % will hold [omega, p1_value] pairs

for om = omega_range

    x = x_star + [0.05; -0.05; 0.02; -0.02];
    diverged = false;

    for t = 1:N_transient
        x = F(x(1), x(2), x(3), x(4), 0, om, om, om, om);
        if any(abs(x) > 1e6) || any(isnan(x))
            diverged = true;
            break
        end
    end

    if diverged
        continue   % skip this omega, too unstable to plot meaningfully
    end

    for t = 1:N_keep
        x = F(x(1), x(2), x(3), x(4), 0, om, om, om, om);
        if any(abs(x) > 1e6) || any(isnan(x))
            break
        end
        bifurcation_p1 = [bifurcation_p1; om, x(1)];
    end
end

figure;
plot(bifurcation_p1(:,1), bifurcation_p1(:,2), '.', 'MarkerSize', 2)
xlabel('\omega (speed parameter)')
ylabel('p1 (long-run values)')
title('Bifurcation diagram: p1 vs adjustment speed')