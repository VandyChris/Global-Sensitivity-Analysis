function y = func_beam(X)
P = X(:, 1);
E = X(:, 2);
v = X(:, 3);
b = X(:, 4);
h = X(:, 5);
L = X(:, 6);

    I = b.*h.^3./12;
    y = P./6./E./I.*((4+5.*v).*h.^2.*L/4+2*L.^3);

end