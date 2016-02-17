function output = beam_RBD( input, Mu, Std )

P = norminv(input(:, 1), Mu(1), Std(1));
E = norminv(input(:, 2), Mu(2), Std(2));
v = norminv(input(:, 3), Mu(3), Std(3));
b = norminv(input(:, 4), Mu(4), Std(4));
h = norminv(input(:, 5), Mu(5), Std(5));
L = norminv(input(:, 6), Mu(6), Std(6));

output = func_beam([P, E, v, b, h, L]);


end

