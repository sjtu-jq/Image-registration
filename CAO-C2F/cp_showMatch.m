
function cp_showMatch(I1,I2,loc1,loc2,correctPos,tname)
cols = size(I1,2);
if size(I1,1) < size(I2,1)
    I1_p = [I1; zeros(size(I2,1)-size(I1,1),size(I1,2),size(I1,3))];
    im3 = [I1_p I2];
elseif size(I1,1) > size(I2,1)
    zeros(size(I1,1)-size(I2,1),size(I2,2));
    I2_p = [I2; zeros(size(I1,1)-size(I2,1),size(I2,2),size(I1,3))];
    im3 = [I1 I2_p];
else
    im3 = [I1 I2];
end

figure,imshow(im3,[])
title(tname);
hold on
if ~isempty(correctPos)
    line([loc1(correctPos,1) loc2(correctPos,1)+cols]',[loc1(correctPos,2) loc2(correctPos,2)]', 'Color', 'y','LineWidth',1.2);
    plot(loc1(correctPos,1),loc1(correctPos,2),'r+');
    plot(loc2(correctPos,1)+cols,loc2(correctPos,2),'go');
    loc1(correctPos,:) = [];
    loc2(correctPos,:) = [];
    line([loc1(:,1) loc2(:,1)+cols]',[loc1(:,2) loc2(:,2)]', 'Color', 'r','LineWidth',1.2);
    plot(loc1(:,1),loc1(:,2),'r+');
    plot(loc2(:,1)+cols,loc2(:,2),'ro');
else
    line([loc1(:,1) loc2(:,1)+cols]',[loc1(:,2) loc2(:,2)]', 'Color', 'y','LineWidth',1.2);
    plot(loc1(:,1),loc1(:,2),'r+');
    plot(loc2(:,1)+cols,loc2(:,2),'go');
end
