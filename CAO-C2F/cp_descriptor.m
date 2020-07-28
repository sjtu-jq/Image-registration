
function descriptors = cp_descriptor(img,keypoints)
NBP = 4;
NBO = 8;
key_num = size(keypoints, 2);
descriptors = zeros([NBP*NBP*NBO key_num]);
%======= compute the image gradient vertors =====%
[M, N] = size(img); % M is the height of image, N is the width of image
du_filter = [-1 0 1];
dv_filter = du_filter';
duv_filter = diag([-1 0 1]);
dvu_filter = flip(duv_filter,2);
gradient_u = imfilter(img, du_filter);
gradient_v = imfilter(img, dv_filter);
gradient_uv = imfilter(img,duv_filter);
gradient_vu = imfilter(img,dvu_filter);
gradient_x = 1.414*gradient_u + gradient_uv - gradient_vu;
gradient_y = 1.414*gradient_v + gradient_uv + gradient_vu;
% dx_filter = [-1 0 1];
% dy_filter = dx_filter';
% gradient_x = imfilter(img, dx_filter);
% gradient_y = imfilter(img, dy_filter);
magnitudes = sqrt( gradient_x.^2 + gradient_y.^2);
%====== compute angle of gradient =======%
angles = zeros(size(img));
for i = 1:M
    for j = 1:N
        if gradient_x(i,j) <= 0
            angles(i,j) = atan(gradient_y(i,j)/gradient_x(i,j)) + pi;
        elseif gradient_x(i,j) > 0 && gradient_y(i,j) >= 0
            angles(i,j) = atan(gradient_y(i,j)/gradient_x(i,j));
        elseif gradient_x(i,j) > 0 && gradient_y(i,j) < 0
            angles(i,j) = atan(gradient_y(i,j)/gradient_x(i,j)) + 2*pi;
        end
    end
end
x = keypoints(1,:); %u
y = keypoints(2,:); %v
s = keypoints(3,:); %scale
theta = keypoints(4,:);
sintheta = sin(theta);
costheta = cos(theta);

[xx,yy] = meshgrid(-NBP/2 : NBP/2);
wincoef = exp( -(xx.^2 + yy.^2)/NBP^2 * 2);
%========= compute descriptors ==========%
parfor p = 1: key_num
    magnitude = magnitudes;
    sp = s(p);
    xp= x(p); % u == x
    yp= y(p); % v == y
    sinth0 = sintheta(p) ;
    costh0 = costheta(p) ;
    W = sp; % scale
    ss = W/2; % scale / 2

    descriptor = zeros(NBP, NBP, NBO); % 4 * 4 * 8
    
    pp = magnitudes(max(-W, 1-yp)+yp : min(+W, M -yp)+yp, max(-W, 1-xp)+xp: min(W, N- xp)+xp);
    pp = pp./max(pp(:));
    xx = sort(pp(:));
    xind1 = round((1 - 1/5) * size(xx,1)); 
    xind2 = round((1 - 2/5) * size(xx,1));
    xind3 = round((1 - 3/5) * size(xx,1));
    xind4 = round((1 - 4/5) * size(xx,1));
    pp = (pp>=xx(xind1)) + (pp<xx(xind1) & (pp>=xx(xind2))) * 0.75 + (pp<xx(xind2) & (pp>=xx(xind3))) * 0.5...
            +(pp<xx(xind3) & (pp>=xx(xind4))) * 0.25;
    magnitude(max(-W, 1-yp) +yp: min(+W, M -yp)+yp, max(-W, 1-xp)+xp: min(W, N- xp)+xp) = pp;

    [dx,dy] = meshgrid(max(-W, 1-xp): min(W, N - xp),max(-W, 1-yp) : min(+W, M - yp));
    nx = (costh0 * dx + sinth0 * dy ) / ss;
    ny = (-sinth0 * dx + costh0 * dy ) / ss;
    for kk = 1:numel(dx)
        mag = magnitude(yp + dy(kk), xp + dx(kk)); 
        angle = angles(yp + dy(kk), xp + dx(kk)) ;  
        angle = mod(angle - theta(p), pi); 
% Cubic interpolation           
        nt = NBO * angle / pi ;
        binx = floor( nx(kk) - 0.5 ) ;
        biny = floor( ny(kk) - 0.5 ) ;
        bint = floor( nt );
        rbinx = nx(kk) - (binx+0.5) ;
        rbiny = ny(kk) - (biny+0.5) ;
        rbint = nt - bint ;
        for dbinx = 0:1
           for dbiny = 0:1
               for dbint = 0:1
                    if  binx+dbinx >= -(NBP/2) && ...
                        binx+dbinx <   (NBP/2) && ...
                        biny+dbiny >= -(NBP/2) && ...
                        biny+dbiny <   (NBP/2) &&  ~isnan(bint)

                          weight = wincoef(binx+dbinx + NBP/2 + 1, biny+dbiny + NBP/2+ 1)...
                              * mag * abs(1 - dbinx - rbinx) * abs(1 - dbiny - rbiny) ...
                              * abs(1 - dbint - rbint) ;

                          descriptor(binx+dbinx + NBP/2 + 1, biny+dbiny + NBP/2+ 1, mod((bint+dbint),NBO)+1) = ...
                              descriptor(binx+dbinx + NBP/2+ 1, biny+dbiny + NBP/2+ 1, mod((bint+dbint),NBO)+1 ) +  weight ;
                    end
               end
           end
        end
    end
    descriptor = reshape(descriptor, 1, NBP * NBP * NBO);
    descriptor = descriptor ./ norm(descriptor);
    
    descriptors(:,p) = descriptor'; % arrange based on column 
end
