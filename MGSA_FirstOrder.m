function index = MGSA_FirstOrder( input, output, ndomain, varargin)
% input: nsample * nd matrix, where nsample is the number of sample, and nd
% is the input dimension
% output: nsample * 1 array
% ndomain: number of sub-domain the to divide a single input, sqrt(nsample)
% recommended and nsample larger than 50 recommended.
% varargin: incudes two arg. arg 1 is a vector of dimension number of
% inputs for which cdf fucntions are available. For example, [2, 3] means
% input 2 and input 3 have cdf functions. arg 2 is a cell of corresponding cdf functions

% This algorithm is proposed by me. Please cite the following paper if you use this code.
% Li, Chenzhao, and Sankaran Mahadevan. "An efficient modularized sample-based method to estimate the first-order Sobol’index." Reliability Engineering & System Safety (2016).


[nsample, nd] = size(input);

%% convert the samples into cdf domains
U = linspace(0, 1, ndomain+1);
U = U';

cdf_input = zeros(nsample, nd);

cdf_values = transpose(linspace(1/nsample, 1, nsample));
if nargin == 3 % cdf function provided
    varargin = {[], {}};
end

j = 1;
for i = 1:nd
    if ismember(i, varargin{1})
        cdf_input(:, i) = varargin{2}{j}(input(:, i));
        j = j+1;
    else
        [~, IX] = sort(input(:, i));
        [~, IX2] = sort(IX);
        cdf_input(:, i) = cdf_values(IX2);
    end
end

%% compute indices
VY = var(output);
VarY_local = zeros(ndomain, nd);
for i = 1:nd
    cdf_input_i = cdf_input(:, i);
    output_i = output;
    U_i = U;
    for j = 1:ndomain
        sub = cdf_input_i<U_i(j+1);
        VarY_local(j, i) = var(output_i(sub));
        inverse_sub = ~sub;
        cdf_input_i = cdf_input_i(inverse_sub);
        output_i = output_i(inverse_sub);
    end
end
index = 1-mean(VarY_local)/VY;
end

