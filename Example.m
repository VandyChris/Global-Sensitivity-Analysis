clc;
clear;
close all
rng(1);
%% Assume uncorrelated inputs
% the samples
Mu = [2.5,200,0.225,1.0,3,3.5]; % mean value
Std = Mu*0.1; % standard deviation
Sigma = diag(Std.^2); % standard deviation;
N = 2e5;
n = 2e3;
input_sample = lhsnorm(Mu, Sigma, N);

% the double loop method
First_order_double_loop   = GSA_FirstOrder(input_sample, @func_beam, n);
Total_effects_double_loop = GSA_TotalEffect(input_sample, @func_beam, n);

% Saltelli's single loop method
A = input_sample;
B = lhsnorm(Mu, Sigma, N);
First_order_Saltelli = Sen_FirstOrder_Saltelli(A, B, @func_beam);
Total_effects_Saltelli = Sen_TotalEffect_Saltelli(A, B, @func_beam);

% Sobol' method
First_order_Sobol = Sen_FirstOrder_Sobol_07(A, B, @func_beam);

% improved FAST based on RBD
k = size(Mu, 2);
First_order_FAST = Sen_FirstOrder_RBD(N, k, @beam_RBD, Mu, Std);

%% plot
figure(1)
plot(First_order_double_loop, 'r+'); hold on;
plot(First_order_Saltelli, 'bs'); hold on;
plot(First_order_Sobol, 'k^'); hold on;
plot(First_order_FAST, 'mv');
xlim([0, 7]);
set(gca, 'XTickLabel',{'P','E','v','b','H','L'},'XTick',[1 2 3 4 5 6]);
xlabel('input')
ylabel('Total effects index')
legend({'Double loop', 'Saltelli', 'Sobol', 'FAST'}, 'location', 'northwest')
title('First-order index, uncorrelated input')

figure(2)
plot(Total_effects_double_loop, 'r+'); hold on;
plot(Total_effects_Saltelli, 'bs'); hold on;
xlim([0, 7]);
set(gca, 'XTickLabel',{'P','E','v','b','H','L'},'XTick',[1 2 3 4 5 6]);
xlabel('input')
ylabel('First-order index')
legend({'Double loop', 'Saltelli'}, 'location', 'northwest')
title('Total effects index, uncorrelated input')

%% Assume correlated inputs
corr = [1	0.174	0.451	0.082	-0.134	0.004;
0.174	1	-0.8	0.059	-0.125	-0.082;
0.451	-0.8	1	-0.004	0.033	0.08;
0.082	0.059	-0.004	1	-0.105	-0.4;
-0.134	-0.125	0.033	-0.105	1	0.279;
0.004	-0.082	0.08	-0.4	0.279	1]; % correlation matrix

COV = zeros(k, k);
for i = 1:k
    for j = 1:k
        COV(i, j) = corr(i, j) * Std(i) * Std(j);
    end
end
First_order_double_loop_2 = GSA_FirstOrder_mvn(Mu', COV, @func_beam, 1000);

figure(3)
plot(First_order_double_loop_2, 'bs'); hold on;
xlim([0, 7]);
set(gca, 'XTickLabel',{'P','E','v','b','H','L'},'XTick',[1 2 3 4 5 6]);
xlabel('input')
ylabel('First-order index')
legend({'Double loop'}, 'location', 'northwest')
title('Total effects index, correlated input')
