function Sen_vector = Sen_FirstOrder_RBD( N, k, model, varargin)
% this function compute the first order Sobol index using the imporved FAST
% method [1] based on random balanced design
% N: sample number
% k: dimension of the inputs
% model: the determinisitic function mapping inputs to outputs
% varargin: parameters of the model
% Note: this functio assumes that all inputs are independent and follows
% the standard uniform distribution U(0, 1). If your inputs follow other
% distributions, in your model for each input use the inverse CDF function 
% to map U(0, 1) to your input distributioin, i.e., U(0, 1) -> input -> ouput
% Note: the computational cost of this function is N
% Ref: Tarantola S, Koda M. Random balance designs for the estimation of 
% first order global sensitivity indices 2006
%%
M = 6;
Sen_vector = zeros(1, k);
s0=[-pi:2*pi/N:pi]';
s = zeros(N, k);
for i = 1:k
s0_random=s0(randperm(N)); %Performs a random permutation of the integers from 1 to N
s(:, i) = s0_random;
end

x=0.5+asin(sin(1*s))/pi; % between 0 and 1
x = min(x, 0.9999);
x = max(x, 0.00001);
y=model(x, varargin{:});
for i = 1:k
[~,index]=sort(s(:, i));
yr=y(index);
spectrum=(abs(fft(yr))).^2/N;
V1=2*sum(spectrum(2:M+1)); 
V=sum(spectrum(2:N));
Sen_vector(i)=V1/V;
display([num2str(i), ' out of ', num2str(k), ' indices finished: ', num2str(Sen_vector(1, i))]);
end


end

