function [ Mu_c, Sigma_c , sample_unfixed, sample_total] = mnv_conditional( Mu, Sigma, index_fixed, value, varargin )
% inputs
% Mu: N*1 column vector of mean values
% Sigma: N*N covariance matrix
% index: indicating which variables are fixed. For instance, index_fixed = [1, 0, 0, 1]
% means x1 and x4 are fixed and q=2 and N=4
% value: column vector of the value of the fixed variables
% varargin: sample numbers if you want to sample from the conditional
% multivariate normal distribution

% outpus
% Mu_c: (N-q)*1 column vector of the mean value of conditional normal distribution
% Sigma_c: (N-q)*(N-q) covariance matrix of conditional normal distribution
% sample_unfixed: samples from the conditional normal distribution, each column
% corresponding to an unfixed variable. Sample is empty if varargin is not
% specified
% sample_total: including the fixed variables

%% partition mean values and covariance matrix
index_fixed = logical(index_fixed);
index_unfixed = (index_fixed == 0);
Mu_1 = Mu(index_unfixed); % mean value of unfixed variables
Mu_2 = Mu(index_fixed); % mean values of fixed variables
Sigma_11 = Sigma(index_unfixed, index_unfixed); % covariance matrix of unfixed variables
Sigma_22 = Sigma(index_fixed, index_fixed); % covariance matrix of fixed variables
Sigma_12 = Sigma(index_unfixed, index_fixed);
Sigma_21 = Sigma(index_fixed, index_unfixed);
%% compute the conditional mean and conditional covariance
Mu_c = Mu_1 + (Sigma_12/Sigma_22) * (value-Mu_2);
Sigma_c = Sigma_11 - (Sigma_12/Sigma_22) *Sigma_21;
%% sample from the conditional distribution if required
if nargin>4;
    n = varargin{:};
    sample_unfixed = lhsnorm(Mu_c, Sigma_c, n);
    sample_total = zeros(n, length(Mu));
    sample_total(:, index_unfixed) = sample_unfixed;
    sample_total(:, index_fixed) = repmat(value', [n, 1]);
else
    sample_unfixed = [];
    sample_total = [];
end
end
