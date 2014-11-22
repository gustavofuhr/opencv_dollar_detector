% 
% It does the conversion from world coordinate system to image coordinate
% using the projection matrix
% 
% USAGE
%  im_point = wcs2ics(w_point, P)
%
function im_point = wcs2ics(w_point, P)

w_point = double(w_point);

% test if the im_point is a homogeneous coordinate
if (size(w_point,1) ~= 4)
    error('The input point is not an homogeneous coordinates');
end

im_point      = P*w_point;
im_point(1,:) = im_point(1,:)./im_point(3,:);
im_point(2,:) = im_point(2,:)./im_point(3,:);
im_point      = im_point(1:2,:);
