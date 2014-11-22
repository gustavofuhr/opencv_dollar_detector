clear;
close all;
clc;


I = imread('pets_background.png');
imshow(I); hold on;


P = 1.0e+07 * ... 
   [0.000092758253702  -0.000079086639724  -0.000015319476449 1.248195750930406; ...
   -0.000009234442978   0.000001295124261  -0.000122402299869 0.622898830135874; ...
    0.000000077919286   0.000000055889203  -0.000000028372202 0.003546929854700];


% plot_groundplane_limits(P, size(I));



area = [-30000 20000; -30000 20000];
H    = [P(:,1:2) P(:,4)];
step = 1000;
plot_grid(area, H, step, [])



% for now, the area will be input by the user to determine which region
% will generate candidates
plot_candidates(P, area, 5000, size(I));