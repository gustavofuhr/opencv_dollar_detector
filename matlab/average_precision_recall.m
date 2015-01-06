function [precision, recall] = average_precision_recall(pstep, p, r)

precision = [];
recall = [];
for rx = pstep:pstep:1.0
%     ind = 

    inside_p = p(r < rx & r > rx-pstep);
    if ~isempty(inside_p)
        mean_p = mean(inside_p);
    else
        mean_p = 0;
    end
    recall = [recall rx];
    precision = [precision mean_p];
    
end    

figure;
plot(recall, precision, '-r');

