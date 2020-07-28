function  cp_showResult(I1rgb,I2rgb,I1gray,I2gray,aff,checkerboard)
% CP_SHOWRESULT 显示配准结果，原图，融合图，仿射图
fprintf('\n【2】Images transformating...\n')

if size(I1rgb,3) == 3 & size(I2rgb,3) == 3
    I1_aff = uint8(zeros(size(I2rgb)));
    I1_aff(:,:,1) = uint8(imwarp(I1rgb(:,:,1),aff,'OutputView',imref2d(size(I2gray))));
    I1_aff(:,:,2) = uint8(imwarp(I1rgb(:,:,2),aff,'OutputView',imref2d(size(I2gray))));
    I1_aff(:,:,3) = uint8(imwarp(I1rgb(:,:,3),aff,'OutputView',imref2d(size(I2gray))));
    figure,
    subplot(121), imshow(I1rgb), title('Source IR image');
    subplot(122), imshow(I2rgb), title('Source VI image');
%     subplot(223), imshow(I1_aff), title('IR Transformation');
    figure, imshow(cp_grayFusion(I1_aff,I2rgb),[]), title('Gray-level Everage Fusion');
    figure, imshow(cp_rgbCheckerBoard(I1_aff,I2rgb,checkerboard)), title('RGB Cheakerboard Mosaic Image');
else
    I1_aff = imwarp(I1gray,aff,'OutputView',imref2d(size(I2gray)));
    I1_aff = double(I1_aff);
    I1_aff = double(I1_aff)/max(max(I1_aff));
    I2gray = double(I2gray);
    I2gray = double(I2gray)/max(max(I2gray));
    figure,
    subplot(221), imshow(I1gray), title('Source image');
    subplot(222), imshow(I2gray), title('Source image');
    subplot(223), imshow(I1_aff), title('I1 Transformation');
    subplot(224), imshow(I1_aff+I2gray,[]), title('Gray-level Everage Fusion');
end

fprintf('【2】Transformation completed!\n')

