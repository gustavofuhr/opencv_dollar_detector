function plot_groundplane_limits(P, im_size)

% Pegue a terceira linha de iH, denotada por [h1 h2 h3]

% A equa??o da linha de fuga no ground plane deve ser

H     = [P(:,1:2) P(:,4)];
H_inv = inv(H);

h = H(3, :); 


% h(1)*u + h(2)*v + h(3) = 0
% h(1)*x + h(2)*y + h(3) = 0
% h(2)*y = -h(1)*x - h(3)
% y = (-h(1)/h(2))*x - (h(3)/h(2))

m = -h(1)/h(2);
b = -h(3)/h(2);

for x = 0:1:im_size(1)
    y = m*x + b;
    plot(x, y, 'xr');
    fprintf('x:%f y:%f\n', x, y);
end
