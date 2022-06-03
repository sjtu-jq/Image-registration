% Author: JIANG Qian
% Email : jiang_qqq@qq.com
% Date  : 2020-03-07
% All rights reserved.
function Imosaic = cp_graymosaic(I1,I2,affmat)
% CP_MOSAIC
% input: a pair of source images
% output: a mosaic image based on liniear interplation transformation
[r1,c1] = size(I1); % IR image
[r2,c2] = size(I2); % visibleimage
Imosaic = zeros([r2+2*max(r1,c1) c2+2*max(r1,c1)]);

%%
affinemat = affmat.T;
[v,u] = find(~isnan(Imosaic));
v = v - max(r1,c1); u = u - max(r1,c1);
utvt = [u v ones(length(v),1)] / affinemat;
ut = utvt(:,1)./utvt(:,3);  vt = utvt(:,2)./utvt(:,3);
utu = reshape(ut,[r2+2*max(r1,c1) c2+2*max(r1,c1)]);
vtv = reshape(vt,[r2+2*max(r1,c1) c2+2*max(r1,c1)]);

Iterp = interp2(double(I1), utu, vtv);
[vn,un] = find(~isnan(Iterp));
vmin1 = min(vn); vmax1 = max(vn);
umin1 = min(un); umax1 = max(un);
Imosaic(:) = NaN;
Imosaic(max(r1,c1)+1:max(r1,c1)+r2, max(r1,c1)+1:max(r1,c1)+c2) = I2;
for i=1:size(Imosaic,1)
    for j=1:size(Imosaic,2)
        if ~isnan(Iterp(i,j)) & ~isnan(Imosaic(i,j))
            Imosaic(i,j) =  (Imosaic(i,j)+Iterp(i,j))/2;
        elseif ~isnan(Iterp(i,j)) & isnan(Imosaic(i,j))
            Imosaic(i,j) = Iterp(i,j);
        end
    end
end
validuv= [min(vmin1,max(r1,c1)) min(umin1,max(r1,c1));max(vmax1,max(r1,c1)+r2) max(umax1,max(r1,c1)+c2)];
Imosaic = uint8( Imosaic(validuv(1):validuv(2),validuv(3):validuv(4))) ;

end

