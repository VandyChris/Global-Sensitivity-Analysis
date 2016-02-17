function sen_vector = Sen_FirstOrder_Saltelli( A, B, f, varargin)
% This function compute the first-order Sobol' index by the single loop
% method proposed by Saltelli [1]
% n: number of sample
% k: number of dimension
% A: n*k matrix, quasi-random sample recommended
% B: n*k matrix, quasi-random sample recommended
% you can generate quasi-random sample using Latin Hypercube sampling 
% f: function handle, so y_A = f(A), y_B = f(B)
% Note: the computational cost of this functio is 2*n + n*k
% Ref: [1] Saltelli, A., Ratto, M., Andres, et. al, "Global sensitivity 
% analysis: the primer", 2008, page 164

y_A = f(A, varargin{:});
[N, k] = size(A);
% V_Y = var([y_A; y_B]);
sen_vector = zeros(1, k);
    function S_i = index_i(order)
        f_0 = sum(y_A)/N;
        i = order;
        C_i = B;
        C_i(:, i) = A(:, i);
        y_C_i = f(C_i, varargin{:});
        S_i = (dot(y_A, y_C_i)/N-f_0^2)/(dot(y_A, y_A)/N - f_0^2);
    end

for j = 1:k
    sen_vector(j) = index_i(j);
    display([num2str(j), ' out of ', num2str(k), ' indices finished: ', num2str(sen_vector(j))]);
end
   
end

