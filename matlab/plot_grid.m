

function plot_grid(area, H, step, trans_xy)

if ~isempty(trans_xy)
    area(1,:) = area(1,:) + trans_xy(1);
    area(2,:) = area(2,:) + trans_xy(2);
end


for x = area(1,1):step:area(1,2)    
    pts = [x x; area(2,1) area(2,2); 1 1];
    pts = H*pts;
    pts(1,:) = pts(1,:)./pts(3,:);
    pts(2,:) = pts(2,:)./pts(3,:);
    plot(pts(1,:), pts(2,:), '-', 'Color', [0.4, 0.4, 0.4]);
end

for y = area(2,1):step:area(2,2)    
    pts = [area(1,1) area(1,2); y y; 1 1];
    pts = H*pts;
    pts(1,:) = pts(1,:)./pts(3,:);
    pts(2,:) = pts(2,:)./pts(3,:);
    plot(pts(1,:), pts(2,:), '-', 'Color', [0.4, 0.4, 0.4]);
end
