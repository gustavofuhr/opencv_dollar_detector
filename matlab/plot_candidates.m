function plot_candidates(P, area, step, im_size)

min_im_height = 50;
w_height      = 1.65;
aspect_ratio  = 0.26;

for x = area(1,1):step:area(1,2)
    for y = area(2,1):step:area(2,2)
        % see if the bounding box is inside the frame
        
        w_point = [x, y];
        bbox = wcoord2bbox(w_point, P, w_height, aspect_ratio);
        
        % for now, only if the whole bbox is inside the image
        % i will consider
        if bbox(1,1) > 0 && bbox(1,2) < im_size(1) && ...
                bbox(1,2) > 0 && bbox(2,2) < im_size(2)
            % check min height
            im_height = bbox(1,2) 
        end
                
    end
end
        


function [bbox] = wcoord2bbox(w_point, P, w_height, aspect_ratio)


% first thing is to find the feet and head points
w_feet_point = [w_point; 0.0; 1.0];
w_head_point = [w_point; w_height; 1.0];

i_feet_point = wcs2ics(w_feet_point, P);
i_head_point = wcs2ics(w_head_point, P);

i_height = abs(i_feet_point(2) - i_head_point(2));
i_width  = i_height * aspect_ratio;

% the middle point is the mean of x (fits diagonal poles better)
ix_middle = (i_feet_point(1) + i_feet_point(1))/2.0;
bbox = [(ix_middle - i_width/2.0) (ix_middle + i_width/2.0); ...
                  i_feet_point(2)           i_head_point(2)];
    