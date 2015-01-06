% this script process the results save in the logFile of the ocd detector
% and generates a precision-recall curves.
clear; close all; clc;


in_file = 'in/calibrated_aspectratio.txt';
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
for f = begin_frame:end_frame
    f
    % only for debugging
    options.image_suf = '/Users/gfuhr/phd/datasets/pets/Crowd_PETS09/S2/L1/Time_12-34/View_001/frame_';
    options.d_mask    = 4;
    options.file_ext  = 'jpg';  
    im_frame = get_frame(options, f+1);
    
    gt_bb = gt{f};
    gt_bb = gt_bb(:, 2:end);
    
    % use different thresholds to sample the curve of precision recall
    fine_th = -1.01:0.05:1;
    for th = fine_th
        bb = filter_detections(all_columns, th, f, resize_factor);
        [precision, recall] = precision_recall_score(gt_bb, bb);
        ps = [ps precision];
        rs = [rs recall];
        
%         imshow(im_frame); hold on;
%         showboxes(im_frame, gt_bb, [0,0,1]);
%         showboxes(im_frame, new_det, [1,0,0]);
    end
    
    coarser_th = -1:2:100;
    for th = coarser_th
        bb = filter_detections(all_columns, th, f, resize_factor);
        [precision, recall] = precision_recall_score(gt_bb, bb);
        ps = [ps precision];
        rs = [rs recall];
        
%         imshow(im_frame); hold on;
%         showboxes(im_frame, gt_bb, [0,0,1]);
%         showboxes(im_frame, new_det, [1,0,0]);
    end
    
    
    
%     pause;
    
end
[fprecision, frecall] = average_precision_recall(0.05, ps, rs);




