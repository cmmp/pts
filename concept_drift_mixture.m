rng default

X = zeros(2000, 2);

% plotmatrix(Z)

t = 1:2000;
w1 = 2 .^ (-0.001 .* t);
w2 = 1 - w1;

for i = 1:2000
   X(i,:) = w1(i) * (5 + randn(1, 2)) + w2(i) * (50 + randn(1,2));
end

slice = 1:100;
scatter(X(slice,1), X(slice,2));