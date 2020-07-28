function grayfusion = cp_grayFusion(I1rgb,I2rgb)
% CP_GRAYFUSION obtain fusion results by everage fusion of two gray images
I1gray= double(rgb2gray(I1rgb));
I1gray = double(I1gray)/max(max(I1gray));
I2gray = double(rgb2gray(I2rgb));
I2gray = double(I2gray)/max(max(I2gray));
grayfusion = I1gray + I2gray;
end

