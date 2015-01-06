clear;
close all;
clc;


load('precision_recall.mat');

mean_p = [];
for rx = 0.1:0.05:1.0
%     ind = 
    mean_p = [mean_p mean(p(r < rx & r > rx-0.05))];
end    

plot([0.1:0.05:1.0], mean_p, '-r');
axis([0, 1.0, 0.0 1.0])