function Icombine = cp_combineImages(numofpart);
%CP_COMBINEIMAGES input images of some parts of the equipment scene, and
%output one complete image of whole equipment
if nargin==0
    numofpart = input('please input the num of images to be combined£º\n');
end
[Filename1,pathname] = uigetfile('*.*','pick 1(st) image');
Icombine = imread([pathname Filename1]);
if size(Icombine,3)==3
    for i = 2:numofpart
    %     [Filename1,pathname] = uigetfile([pathname '*' Filename1(end-3:end)],'pick next image');
        Ipart = imread([pathname Filename1(1:end-5) num2str(i) Filename1(end-3:end)]) ;
        [x,y] = find(Ipart(:,:,1)~=0 | Ipart(:,:,2)~=0 | Ipart(:,:,3)~=0);
        for j=1:length(x)
            if sum(Icombine(x(j),y(j),:))==0
                Icombine(x(j),y(j),:) =Ipart(x(j),y(j),:);
            elseif sum(Icombine(x(j),y(j),:))~=0 && (sum(Ipart(x(j),y(j),:)) > sum(Icombine(x(j),y(j),:)+20))
                Icombine(x(j),y(j),:) = Ipart(x(j),y(j),:);
            end
        end
    %     Icombine(:) = max(Icombine(:), Ipart(:));
    end
else
    for i = 2:numofpart
    %     [Filename1,pathname] = uigetfile([pathname '*' Filename1(end-3:end)],'pick next image');
        Ipart = imread([pathname Filename1(1:end-5) num2str(i) Filename1(end-3:end)]) ;
        [x,y] = find(Ipart(:,:)>0);
        for j=1:length(x)
            if Icombine(x(j),y(j))==0
                Icombine(x(j),y(j)) =Ipart(x(j),y(j));
            elseif Icombine(x(j),y(j))>0 %&&  ( Ipart(x(j),y(j))> (Icombine(x(j),y(j))+15))
                Icombine(x(j),y(j)) = max(Icombine(x(j),y(j)),Ipart(x(j),y(j)));
            end
        end
    end
end
end
