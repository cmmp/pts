rng default

X = zeros(2000, 2);

% plotmatrix(Z)

t = 1:2000;
w1 = 2 .^ (-0.00002 .* t);
w2 = 1 - w1;

for i = 1:2000
    r = rand();
    if r <= w1(i)
        X(i,:) = 5 + randn(1, 2);
    else
        X(i,:) = 50 + randn(1,2);
    end
end

slice = 1000:1500;
scatter(X(slice,1), X(slice,2));





