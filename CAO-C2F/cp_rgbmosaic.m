function Imosaic = cp_rgbmosaic(I1,I2,affmat)
% CP_MOSAIC
% input: a pair of source rgb images
% output: a mosaic image based on liniear interplation transformation
Imosaic(:,:,1) = cp_graymosaic(I1(:,:,1),I2(:,:,1),affmat);
Imosaic(:,:,2) = cp_graymosaic(I1(:,:,2),I2(:,:,2),affmat);
Imosaic(:,:,3) = cp_graymosaic(I1(:,:,3),I2(:,:,3),affmat);

end

