% this script process the results save in the logFile of the ocd detector
% and generates a precision-recall curves.
clear; close all; clc;


in_file = 'in/calibrated_biggerstd.txt';
ground_truth_detections = 'in/pets_gt_MOT_all.mat';
resize_factor = 1.5; % this is the resize factor used in ocd


% begin the miracle
fid = fopen(in_file);
[all_columns] = textscan(fid, '%f %f %f %f %f %f');

load(ground_truth_detections);

begin_frame = min(all_columns{1});
end_frame   = max(all_columns{1});


ps = [];
rs = [];
        
% % use different thresholds to sample the curve of precision recall
% fine_th = -1.01:0.02:1;
% for th = fine_th
%     fprintf('Processing %d of %d threshold\n', find(fine_th==th), length(fine_th)); tic;
%     bb = filter_detections_th(all_columns, th, resize_factor);
%     
%     [precision, recall] = precision_recall_score_allframes(gt, bb);
%     ps = [ps precision];
%     rs = [rs recall];
%     toc;
% end

coarser_th = 1:2:150;
for th = coarser_th
    fprintf('Processing %d of %d threshold\n', find(coarser_th==th), length(coarser_th)); tic;
    bb = filter_detections_th(all_columns, th, resize_factor);
    
    [precision, recall] = precision_recall_score_allframes(gt, bb);
    ps = [ps precision];
    rs = [rs recall];
    toc;
end

    
[fprecision, frecall] = average_precision_recall(0.01, ps, rs);

new_precision = fprecision(fprecision ~= 0.0);
new_recall = frecall(fprecision ~= 0.0);
plot(new_recall, new_precision, '-b');




