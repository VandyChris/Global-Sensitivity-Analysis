function sen_vector = Sen_FirstOrder_Sobol_07( x, z, f, varargin)
% This function compute the first-order Sobol' index based on the algorithm
% proposed by Sobol' [1]. This function shows higher accuracy than
% Sen_First_Saltelli for small indices
% n: number of sample
% k: number of dimension
% x: n*k matrix, quasi-random sample recommended
% z: n*k matrix, quasi-random sample recommended
% you can generate quasi-random sample using 
% f: function handle, so y_x = f(x), y_z = f(z)
% [1] Sobol' IM, Myshetskaya EE. Monte Carlo estimators for small sensitivity 
% indices. Monte Carlo Methods Appl 2008;13:455–65.

y_x = f(x, varargin{:});
y_z = f(z, varargin{:});
[~, k] = size(x);
Vy = var([y_x; y_z]);
c = mean([y_x; y_z]);
sen_vector = zeros(1, k);

for i = 1:k
    xz = z;
    xz(:, i) = x(:, i);
    y_xz = f(xz, varargin{:});
    V_i = mean((y_x - c).*(y_xz - y_z));
    sen_vector(i) = V_i / Vy;
    display([num2str(i), ' out of ', num2str(k), ' indices finished: ', num2str(sen_vector(i))]);
end

end

