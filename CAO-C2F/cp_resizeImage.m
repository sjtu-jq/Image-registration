
function [Ir1, Ir2, scale] = cp_resizeImage(I1,I2,height)
% ������ͼƬ�ֱ��ʵ���Ϊ����С�߶�Ϊ��׼��ͼ�����Ϊԭʼ����
disp('��Images are Resized!��')
scale =[ size(I1,1) / height ; size(I2,1) / height ];
Ir1 = imresize(I1, height/size(I1,1));
Ir2 = imresize(I2, height/size(I2,1));
end
% h = fspecial('gaussian',3,1.5);
% Ir1 = conv2(Ir1,h,'same');
% Ir2 = conv2(Ir2,h,'same');