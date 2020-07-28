
function [Ir1, Ir2, scale] = cp_resizeImage(I1,I2,height)
% 将两幅图片分辨率调整为以最小高度为基准，图像比例为原始比例
disp('【Images are Resized!】')
scale =[ size(I1,1) / height ; size(I2,1) / height ];
Ir1 = imresize(I1, height/size(I1,1));
Ir2 = imresize(I2, height/size(I2,1));
end
% h = fspecial('gaussian',3,1.5);
% Ir1 = conv2(Ir1,h,'same');
% Ir2 = conv2(Ir2,h,'same');