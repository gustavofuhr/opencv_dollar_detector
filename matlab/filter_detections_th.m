function [detections] = filter_detections_th(all_data, threshold, resize_factor)

min_frame = min(all_data{1});
max_frame = max(all_data{1});


detections = {}; 
for f = min_frame:max_frame
    
    det_frame_th{1} = all_data{2}(all_data{1} == f & all_data{6} > threshold);
    det_frame_th{2} = all_data{3}(all_data{1} == f & all_data{6} > threshold);
    det_frame_th{3} = all_data{4}(all_data{1} == f & all_data{6} > threshold);
    det_frame_th{4} = all_data{5}(all_data{1} == f & all_data{6} > threshold);


    det = [det_frame_th{1}/resize_factor det_frame_th{2}/resize_factor ...
           det_frame_th{3}/resize_factor det_frame_th{4}/resize_factor];
       
    if ~isempty(det)
        det(:,3) = det(:,3) + det(:,2);
        det(:,4) = det(:,4) + det(:,1);
        new_det = det;
        new_det(:,3) = new_det(:,4);
        new_det(:,4) = det(:,3);
        
        detections{f} = new_det;
    else
        detections{f} = [];
    end
    
end