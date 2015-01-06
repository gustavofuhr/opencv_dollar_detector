function [precision, recall] = precision_recall_score(gt_bb, bb);

CompFunctionName = 'VOCscore';
debug = false;

%% Init structures
truepos = 0;
falseneg = 0;
falsepos = 0;

%% distances
distances = 0.;
% gt = 0;

%% Tmp Counter
%trueposTmp = 0;
%falsenegTmp = 0;
%falseposTmp = 0;
Ass = [];
Cost = [];

%% Compute current mapping procedure.
score = [];


%% Getting the distance matrix Tracks vs Annotations
currentAllLabel = [];
for b=1:size(gt_bb,1);
    for l=1:size(bb,1)
        if strcmp(CompFunctionName, 'DistFeet3D')
            distance = mydistance(gt_bb(b,:),bb(l, :), CompFunctionName, H);
        else
            distance = mydistance(gt_bb(b,:),bb(l, :), CompFunctionName);
        end
        score(b,l) = distance;
    end
end

%% From distance matrix get association matrix
if strcmp(CompFunctionName, 'DistFeet3D')
    dist = 500;
    Ass = GreedyAssociation(score,dist,true);
else
    dist = 0.5;
    Ass = GreedyAssociation(score,dist,false);
end;


%% Compute current mapping (between tracks hyp. and annotations).
% Note in mapping there is: [idAnnotation, idTrackerHyp.] where
% idTrackerHyop is the index inside idxTracks().
mapping = [];
for r=1:length(Ass(:))
    if Ass(r) == 1
        [b, l] = ind2sub(size(Ass),r);
        obj = b;
        tt = l;
        mapping = [mapping; obj tt];
    end
end


if length(mapping) > 0
    for o=1:length(mapping(:,1))
        %% count as TP and evaluate the MOTP.
        truepos = truepos + 1;
        %h = find(idxTracks == mapping(o,2));
        %idxo= find(indexObj(:,h) == mapping(o,1));
        %distances = distances + score(idxo,h);
        %trueposTmp = trueposTmp + 1;
    end
end

%% Check false negative (unmapped annotated obj. up to a threshold).
for r=1:size(Ass,1)
    if sum(Ass(r,:)) == 0
        if debug, fprintf('False negative detected for target %d!\n', r); end;
        falseneg = falseneg + 1;
        %falsenegTmp = falsenegTmp + 1;
    end
end
if isempty(Ass)
    for r=1:size(gt_bb,1)
        if debug, fprintf('False negative detected!\n'); end;
        falseneg = falseneg + 1;
        %falsenegTmp = falsenegTmp + 1;
    end
end

%% Check false positive (unmapped tracker hyp. up to a threshold).
for c=1:size(Ass,2)
    if sum(Ass(:,c)) == 0
        if debug, fprintf('False positive detected for target %d!\n', c); end;
        falsepos = falsepos +1;
        %falseposTmp = falseposTmp + 1;
    end
end

if (truepos + falsepos) ~= 0
    precision = double(truepos) / (truepos + falsepos);
else
    precision = 0.0;
end

if (truepos + falseneg) ~= 0
    recall = double(truepos) / (truepos + falseneg);
else
    recall = 0.0;
end


end

function dist = mydistance(bboxesDetect, target, typeComp, H)

if(strcmp('VOCscore',typeComp))
    xtlA = bboxesDetect(1,1);
    ytlA = bboxesDetect(1,2);
    woA = bboxesDetect(1,3)-bboxesDetect(1,1);
    hoA = bboxesDetect(1,4)-bboxesDetect(1,2);
    
    
    xtlT =  target(1);
    ytlT =  target(2);
    woT =  target(3)-target(1);
    hoT =  target(4)-target(2);
    
    intersection = rectint([ xtlT ytlT woT hoT ],[ xtlA ytlA woA hoA ]);
    union = (woT*hoT) + (woA*hoA) - intersection;
    if(union == 0)
        dist=0;
    else
        dist = intersection/union;
    end
elseif (strcmp('DistFeet3D',typeComp))
    % if this is chosen, the user should pass also an homography H
    gt_tl_x = bboxesDetect(1,2);
    gt_br_x = bboxesDetect(1,4);
    gt_br_y = bboxesDetect(1,5);
    
    gt_feet = inv(H)*[(gt_tl_x + gt_br_x)/2; gt_br_y ; 1];
    gt_feet(1) = gt_feet(1)/gt_feet(3); gt_feet(2) = gt_feet(2)/gt_feet(3);
    gt_feet = gt_feet(1:2);
    
    t_tl_x =  target.bbox(1);
    t_br_x = t_tl_x + target.bbox(3);
    t_br_y = target.bbox(2) + target.bbox(4);
    
    t_feet = inv(H)*[(t_tl_x + t_br_x)/2; t_br_y ; 1];
    t_feet(1) = t_feet(1)/t_feet(3); t_feet(2) = t_feet(2)/t_feet(3);
    t_feet = t_feet(1:2);
    
    dist = norm(gt_feet - t_feet);
end
return
end