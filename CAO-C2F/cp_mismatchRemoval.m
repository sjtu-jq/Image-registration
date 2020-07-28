function correctIndex = cp_mismatchRemoval(p1,p2,edge1,edge2,maxErr)
% mismatch removal algorithm
% Input: cp_mismatchRemoval(p1,p2,edge1,edge2,maxErr,iteration)
% p1: coarse match of image 1
% p2: coarse match of image 2
if nargin == 5
    iteration = 1;
    zoomflag = 1;
end
if length(p1) == 3
    correctIndex = [1 2 3];
    regis2sub = p2;
    return;
end
if length(p1) == 2
    correctIndex = [1 2];
    return;
end
correctIndex = [];
minArea1 = size(edge1,1) * size(edge1,2) / 80;
minArea2 = size(edge2,1) * size(edge2,2) / 80;

len = size(p1,1);
eor = [];
%% RANSAC algorithm
% if  ~zoomflag
%     for i = 1:len-2
%         for j = i+1:len-1
%             for k = j+1:len
%                 pp1 = [p1(i,1:2);p1(j,1:2);p1(k,1:2)];
%                 pp2 = [p2(i,1:2);p2(j,1:2);p2(k,1:2)];
%                 A1 = p1(i,:)-p1(j,:);
%                 B1 = p1(i,:)-p1(k,:);
%                 A2 = p2(i,:)-p2(j,:);
%                 B2 = p2(i,:)-p2(k,:);
%                 if 4*abs( A1(1)*B1(2)-A1(2)*B1(1) ) < minArea1 || 4*abs( A2(1)*B2(2)-A2(2)*B2(1) ) < minArea2
%                     continue;
%                 end
%                 side1A = norm(A1); side1B = norm(B1);
%                 side2A = norm(A2); side2B = norm(B2);
%                 A2A = side1A / side1B;
%                 B2B = side2A / side2B;
%                 if  A2A/B2B < 0.92 || A2A/B2B > 1.08
%                     continue;
%                 end
%                 ERR = getEdgeOverlappingRatio(pp1,pp2,p1(:,1:2),p2(:,1:2),i,j,k,maxErr);
%                 eor = [eor;i j k ERR]; 
%             end
%         end
%     end
% else
    for n = 1:min(len*(len-1)*(len-2)/6,500) 
        ijk = randperm(len,3);
        ir = ijk(1);
        jr = ijk(2);
        kr = ijk(3);
        pp1 = [p1(ir,1:2); p1(jr,1:2); p1(kr,1:2)];
        pp2 = [p2(ir,1:2); p2(jr,1:2); p2(kr,1:2)];
        A1 = p1(ir,:)-p1(jr,:);
        B1 = p1(ir,:)-p1(kr,:);
        A2 = p2(ir,:)-p2(jr,:);
        B2 = p2(ir,:)-p2(kr,:);
        if abs( A1(1)*B1(2)-A1(2)*B1(1) ) < minArea1 || abs( A2(1)*B2(2)-A2(2)*B2(1) ) <minArea2
            continue;
        end
        ransacERR = getEdgeOverlappingRatio(pp1,pp2,p1(:,1:2),p2(:,1:2),ir,jr,kr,maxErr);
        eor = [eor; ir jr kr ransacERR];
    end
% end
%% 
if size(eor,1) < 3
    error('No sufficent matches! ( less than 3 matches obtained )');
end
[~,ind] = max(eor(:,5));
base1 = p1(eor(ind,1:3),1:2);
base2 = p2(eor(ind,1:3),1:2);
affmat0 = fitgeotrans(base1,base2,'affine');
correctIndex = eor(ind,1:3)';
for i = 1:len
    if i == eor(ind,1) | i == eor(ind,2) | i == eor(ind,3)
        continue;
    end
    pp1_aff = [p1(i,1:2) 1] * affmat0.T;
    pp2_to_pp1aff = abs(p2(i,1:2) - pp1_aff(1:2));
    if (pp2_to_pp1aff(1) < maxErr / 1.5 & pp2_to_pp1aff(2) < 1.5 * maxErr) | ...
           (pp2_to_pp1aff(1) < 1.5 * maxErr & pp2_to_pp1aff(2) < maxErr / 1.5)
        correctIndex = [correctIndex; i];
    end
end
end
%% RANSAC algorithm
function RMSE = getEdgeOverlappingRatio(pp1, pp2, p1, p2, a, b, c, maxErr)
affmat = fitgeotrans(pp1,pp2,'affine');
p1_aff = [p1 ones(length(p1),1)] * affmat.T;

p2_to_p1aff = (p2 - p1_aff(:,1:2)).^2; % Distance between origin points and affined points
RMSE(1,1) = sqrt(sum(p2_to_p1aff(:)) / length(p2));

pixelDistance = p2_to_p1aff(:,1) + p2_to_p1aff(:,2);
pixelDistance([a b c]) = 2*maxErr^2; % ignore origin three points
[u,~] = find(pixelDistance <  maxErr^2); % mount of internal points of fitting result
[minerr,ind0] = min(pixelDistance );
RMSE(1,2:4) = [length(u) minerr ind0];
end
