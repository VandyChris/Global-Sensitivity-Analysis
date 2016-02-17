function index= GSA_FirstOrder_mvn(Mu, Sigma,func, n, varargin)
% this function compute the first-order Sobol' index for joint Gaussian
% inputs using the double-loop Monte Carlo method
% Mu: column vector of mean values
% Sigma: covariance matrix
% func: a function handle to generate output, which can transfer N*D input
% to N*1 outputs
% n: the number of input samples you'd like to use.
% vargrgin: the argument of func besides input
% Note: high-order Sobol' indices are meaningless for correlated inputs,
% but first-order index is still an informative sensitivity indicator
%%
input=mvnrnd(Mu, Sigma, n);
var_output=var(func(input, varargin{:})); %the variance of the output
D = size(input, 2);
%%
index=zeros(1,D);
for i=1:D;
    output_xi=zeros(n,n);
    index_fixed = zeros(1, D);
    index_fixed(i) = 1;
    for j=1:n;
        func_j = func;
        varargin_j = varargin;
        [~, ~, ~, input_temp] = mnv_conditional(Mu, Sigma, index_fixed, input(j, i), n);
        output_xi(:,j)=func_j(input_temp, varargin_j{:});
    end
    index(1,i)=var(mean(output_xi))/var_output;
    display([num2str(i), ' out of ', num2str(D), ' indices finished: ', num2str(index(1, i))]);
end


end