% CFO-C2F registration for infrared and visible images of power equipment
% Author：Qian Jiang (in chinese:姜骞)
% First version：2019-04-29
% Current version：2020-02-28
%% section I: Read source images 读取源图片
clear all
set(0,'defaultfigurecolor','w') 
DistortFlag = 0;%input('Is there distortion of infrared image? :\n');
eva=[];
for name=[1 2 3 4 7 8 9 10 11 12]
if name<3
    infrared = name;
    visible = name; 
elseif 3<=name && name<=6
    infrared = name;
    visible = 3;
elseif name==7
    infrared = name;
    visible = 4;
elseif 8<=name && name<=12
    infrared = name;
    visible = name;
elseif 13<=name && name<=15
    infrared = name;
    visible = 13;
end
[I1gray, I2gray, I1rgb, I2rgb, f1, f2, path] = cp_readImage;...
%     (0,'D:\研究文件\CSS-SIFT\Exprimental images in Paper\',['I' num2str(infrared) '.jpg'], ['V' num2str(visible) '.jpg']);
%% section III: Resize images based on the minimum imaclosege height
height = size(I1gray,1);
[I1, I2, scale] = cp_resizeImage(I1gray,I2gray,height);
%% section IV: Registrate iteratively & Coarse matching迭代配准
close all;
clc;
I1_itea = I1;
if name==7 || name==15 || name==14
    iterationNum = 2;
else
    iterationNum = 1;
end
iteration = 0;
Runtime = 0;
maxRMSE = 4*ceil(size(I2,1)/300);
AffineTrans = zeros([3 3 iterationNum]);
while  iteration < iterationNum
    fprintf('\n%d(th) iteration of registration...\n',iteration);
    [P1,P2, Rt,corner12] = cp_registration(I1_itea,I2, 20, maxRMSE,iteration, 1,  0,      6, 1    ,I2gray);
                        % cp_registration(I1,    I2, theta,maxRMSE,iteration,zoom+,zoom-,Lc,showflag,I2gray)
    Runtime = Rt + Runtime
    [I1_itea,affmat] = cp_getAffine(I1_itea,I2,P1,P2); % [v1,u1]==[v2,u2]
    iteration = iteration+1;
    AffineTrans(:,:,iteration) = affmat.T;
end
% Points of I1gray after resize 
P1  = [P1 ones([length(P1) 1])];
[pos_cor1,~] = find(corner12(:,1) == 0);
for iteration = iteration:-1:2
    P1 = P1 / AffineTrans(:,:,iteration-1);
    cor12 = [corner12(1:pos_cor1-1,1:2) ones(pos_cor1-1,1)] / AffineTrans(:,:,iteration-1);
    P1(:,1:2) = P1(:,1:2) ./ P1(:,3);
    P1(:,3) = ones(length(P1),1);
    corner12(1:pos_cor1-1,1:2) = cor12(:,1:2) ./ cor12(:,3);
    corner12(1:pos_cor1-1,3) = ones(pos_cor1-1,1);
end
P1 = P1(:,1:2);
corner12 = corner12(:,1:2);
% Correct matches in the source images
P1(:,2) = size(I1gray,1) / 2 + scale(1) * ( P1(:,2)-size(I1,1)/2);
P1(:,1) = size(I1gray,2) / 2 + scale(1) * ( P1(:,1)-size(I1,2)/2);
corner12(1:pos_cor1-1,2) = size(I1gray,1) / 2 + scale(1) * ( corner12(1:pos_cor1-1,2)-size(I1,1)/2);
corner12(1:pos_cor1-1,1) = size(I1gray,2) / 2 + scale(1) * ( corner12(1:pos_cor1-1,1)-size(I1,2)/2);

P2(:,2) = size(I2gray,1) / 2 + scale(2) * ( P2(:,2)-size(I2,1)/2);
P2(:,1) = size(I2gray,2) / 2 + scale(2) * ( P2(:,1)-size(I2,2)/2);
corner12(pos_cor1+1:end,2) = size(I2gray,1) / 2 + scale(2) * ( corner12(pos_cor1+1:end,2)-size(I2,1)/2);
corner12(pos_cor1+1:end,1) = size(I2gray,2) / 2 + scale(2) * ( corner12(pos_cor1+1:end,1)-size(I2,2)/2);
%% section V: Fine matching
tic;
P3 = cp_subpixelFine(P1,P2); % Fine matching
toc
%% section VI: Evaluation critea of the registration results
if exist([path f1(1:end-4) '.mat'],'file')
    load([path f1(1:end-4) '.mat'])
    [Evaluation1,correctindex1] = cp_Evaluation(P1,P2,refTrans,corner12);
    [Evaluation2,correctindex2] = cp_Evaluation(P1,P3,refTrans,corner12);
else
    [Evaluation,correctindex] = cp_Evaluation(P1);
end

%********* Show visual registration result*******%

[~,affmat] = cp_getAffine(I1gray,I2gray,P1,P3);
Imosaic = cp_graymosaic(I1gray, I2gray, affmat);
figure, subplot(121),imshow(Imosaic);subplot(122),imshow(cp_rgbmosaic(I1rgb,I2rgb,affmat));
% cp_showResult(I1rgb,I2rgb,I1gray,I2gray,affmat,3); % checkborder image
% cp_showMatch(I1rgb,I2rgb,P1,P2,[],'Before Subpixel Fining');
cp_showMatch(I1rgb,I2rgb,P1,P3,[],'After Subpixel Fineing');
eva = [eva; name*ones([1 8]); Evaluation1{2}; Evaluation2{2}];
% imwrite(cp_rgbmosaic(I1rgb,I2rgb,affmat),['D:\研究文件\MinePaper\电网技术\Firuge in Paper\' f1(1:end-4) '_Mosaic.jpg']);
end
eva
%% section VII: Write images that have already been registrated
% imwrite(Imosaic,[path f1(1:end-4) '_Mosaic.jpg']);
% Iaffine = imwarp(I1rgb,affmat,'Outputview',imref2d(size(I2gray)));
% imwrite(Iaffine,[path f1(1:end-4) '_Trans.jpg']);
% Icombine = cp_combineImages();figure,imshow(Icombine);
% if exist([path f1(1:end-4) '.mat'],'file')
%     load([path f1(1:end-4) '.mat'])
%     Projective =  [P1 ones(length(P1),1)] * affmat.T; Projective = Projective(:,1:2) ./ Projective(:,3); 
%     Pref =  [P1 ones(length(P1),1)] * refTrans; Pref = Pref(:,1:2) ./ Pref(:,3); 
%     figure;
%     imshow(I1gray);
%     hold on; plot(P1(:,1),P1(:,2),'g+');hold off;
%     figure;imshow(I2gray);
%     hold on;
%     plot(P3(:,1),P3(:,2),'g+');
%     plot(Pref(:,1),Pref(:,2),'ro');
%     plot(Projective(:,1),Projective(:,2),'y*');
% end
%% Obtain reference transformation matrix manually
% [I1_aff,refaffmatT] = cp_manuallyTrans(I1rgb,I2rgb);
% cp_showResult(I1rgb,I2rgb,I1gray,I2gray,refaffmatT,5);
% refTrans = refaffmatT.T;
% save([path f1(1:end-4) '.mat'],'refTrans')