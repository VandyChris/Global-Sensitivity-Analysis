function index= GSA_TotalEffect(input,func,n, varargin)
% this function compute the total effects Sobol' index using the double loop
% Monte Carlo method
% input: a N*D matrix, where N is the number of samples, and D is the 
% dimension of the input
% func: a function to generate output, which can transfer N*D input to N*1 
% outputs
% n: the number of input samples you'd like to use in each loop. For 
% example, N can be 10000 but only 5000 are used in the analysis if you 
% set n=5000
% vargrgin: the argument of func besides input
% Note: this function assumes uncorrelated inputs; the computational cost
% of this function is n^2*D + n, which is pretty high; more efficient
% algorithm can be found in the same folder
%%

[N,D]=size(input);
if N<n;
    disp('Error: not enough sampling points');
    return;
end

input=input(1:n,:);
var_output=var(func(input, varargin{:})); %the variance of the output


%%
index=zeros(1,D);

for i=1:D;
    input_i = input(:, i);
    output_xi=zeros(n,n);
    parfor j=1:n; 
        func_j = func;
        varargin_j = varargin;
        input_j = input(j, :);
        input_temp=repmat(input_j,[n,1]);
        input_temp(:,i)=input_i;
        output_xi(:,j)=func_j(input_temp, varargin_j{:}); 
    end
    index(1,i)=mean(var(output_xi))/var_output;
    display([num2str(i), ' out of ', num2str(D), ' indices finished: ', num2str(index(1, i))]);
end


end