
function [correctIndex1, correctIndex2, zoom] = cp_match(des1, des2, distRatio)
%%
zoomvote = zeros(3,size(des2,3));
match1 = zeros(size(des1,1),size(des2,3));
match2 = zeros(size(des2,2),size(des2,3));
% For each descriptor in the first image, select its match to second image.
for i = 1 : size(des1,1)
    for j = 1:size(des2,3)
       des2t = des2(:,:,j); 
       dotprods = des1(i,:) * des2t;        % Computes vector of dot products
       [vals1,indx1] = sort(acos(dotprods)); % Take inverse cosine and sort results
       % Check if nearest neighbor has angle less than distRatio times 2nd.
       zoomvote(1,j) = 10; % > 2*pi is OK
       if vals1(1) < distRatio * vals1(2)
          zoomvote(1,j) = vals1(1);
          match1(i,j) = indx1(1); % store the matches of each zoom time
       end
    end
%     zoomvote
    [mintheta,ind] = min(zoomvote(1,:));
    if mintheta ~= 10
        zoomvote(2,ind) = zoomvote(2,ind)+1;
    end
end
disp(zoomvote(2,:));
% [~,zoom] =  max(zoomvote(2,:));
% determine zoom scope 
% des2zoom = des2(:,:,zoom)';
% match = match(:,zoom);
%%
des1t = des1';
for j = 1:size(des2,3) 
    des2zoom = des2(:,:,j)';
    for i = 1 : size(des2zoom,1)
        dotprods = des2zoom(i,:) * des1t;        
       [vals1,indx1] = sort(acos(dotprods));
       if vals1(1) < distRatio * vals1(2) 
          match2(i,j) = indx1(1);
       end
    end
end
%% bilateral match
zoomvote(3,:) = size(des1,1) * ones(1,size(des2,3)); %assume all of match1 are correct 
for j = 1:size(des2,3)
    for i = 1 : size(des1,1)
        if match1(i,j)>0 
            if match2(match1(i,j),j)~=i
                match1(i,j)=0;
                zoomvote(3,j) = zoomvote(3,j) - 1;
            end
        else
            zoomvote(3,j) = zoomvote(3,j) - 1;
        end
    end
end
disp(zoomvote(3,:));
[maxvote, zoom] =  max(zoomvote(3,:));
%% obtain match index
correctIndex1 = [];
correctIndex2 = [];
for i = 1 : size(des1,1)
        if match1(i,zoom)>0
            correctIndex1 = [correctIndex1;i];
            correctIndex2 = [correctIndex2;match1(i,zoom)];
        end
end
