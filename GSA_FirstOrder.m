function index= GSA_FirstOrder(input,func,n, varargin)
% this function conducts the sensitivity analysis using the first-order
% index

%here input is a N*D matrix, where N is the number of sampling points, and
%D is the dimension of the input
%func is a function handle to generate output, which can transfer N*D input
%to N*1 outputs, such as a GP model
%n is the number of input samples you'd like to use. For example, N can be
%10000 but only 5000 are used in the analysis if you set n=5000

% vargrgin: the argument of func besides input

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
    output_xi=zeros(n,n);
    parfor j=1:n;
        func_j = func;
        varargin_j = varargin;
        input_temp=input;
        input_temp(:,i)=ones(n,1)*input(j,i);
        output_xi(:,j)=func_j(input_temp, varargin_j{:});
    end
    index(1,i)=var(mean(output_xi))/var_output;
    display([num2str(i), ' out of ', num2str(D), ' indices finished: ', num2str(index(1, i))]);
end


end