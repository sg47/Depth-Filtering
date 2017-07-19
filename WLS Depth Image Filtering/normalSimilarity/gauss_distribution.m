function f = gauss_distribution(x, mu, sCov)
data(1) = x(1,1,1); data(2) = x(1,1,2); data(3) = x(1,1,3);
s(1) = sCov(1,1); s(2) = sCov(2,2); s(3) = sCov(3,3);
p1 = -.5 * ((data - mu)/s) .^ 2;
p2 = (s * sqrt(2*pi));
f = abs(mean(exp(p1) ./ p2)); 