function undistimg=cp_undistortimg(inputimg,internalmat,distort)
%张正友畸变矫正法：去除红外图片的畸变,输入参数为畸变图像，畸变参数，内参矩阵
%%
row = size(inputimg,1);
col = size(inputimg,2);
color = size(inputimg,3);
%%
inter1 = [828.06946 0 242.69122
          0  830.54110 382.53820
          0 0 1]; 
inter2 = [1064.62006 0 287.45702
          0 1065.38040 397.45537
          0 0 1];
dist1 = [-0.40655   0.11804   0.00323   0.00735];
dist2 = [0.24799   -1.17408   0.00006   0.00107];
if nargin==1
    if size(inputimg,1) == 768 && size(inputimg,2) == 576
        internalmat = inter1;
        distort = dist1;
    elseif size(inputimg,1) == 800 && size(inputimg,2) == 600
        internalmat = inter2;
        distort = dist2;
    else
%         distort = [-4.6e-7 0 0 0];
%         internalmat = [1 0 col/2;0 1 row/2;0 0 1];
        undistimg = inputimg;
        return;
    end
end
% if nargin==1
%     distort = [-4.6e-7 0 0 0];
%     internalmat = [1 0 col/2;0 1 row/2;0 0 1];
% elseif nargin == 2
% %     internalmat = [1 0 col/2;0 1 row/2;0 0 1];
%     distort = [-4.6e-7 0 0 0];
% elseif nargin==3
%     if isempty(internalmat)
%         internalmat = [1 0 col/2;0 1 row/2;0 0 1];
%     end
%     if isempty(distort)
%         distort = [-4.6e-7 0 0 0];
%     end
% end

%%
if color==1
    undistimg = im2double(zeros([row col]));
else
    undistimg = uint8(zeros([row col color]));
end
%%
%畸变参数矩阵
k1 = distort(1); k2 = distort(2);
p1 = distort(3); p2 = distort(4);
fu = internalmat(1,1); fv = internalmat(2,2);
u0 = internalmat(1,3); v0 = internalmat(2,3);
%%
%创建坐标网格
    [v,u]=find(~isnan(undistimg(:,:,1)));
%%
%计算uv到xy的变换值和中心点的值
    xyz = [u v ones(row*col,1)] / internalmat';
    r2 = xyz(:,1).^2 + xyz(:,2).^2;
    x = xyz(:,1); y = xyz(:,2);
%%  
%计算对应的畸变前的坐标
    x = x.*(1+k1*r2 + k2*r2.^2) + 2*p1.*x.*y + p2*(r2 + 2*x.^2);
    y = y.*(1+k1*r2 + k2*r2.^2) + 2*p2.*x.*y + p1*(r2 + 2*y.^2);
    u_d = reshape(x * fu + u0,[row col]);
    v_d = reshape(y * fv + v0,[row col]);
%%
%双线性插值
    if color == 1
        undistimg = interp2(double(inputimg),u_d,v_d);
    elseif color == 3
        R = double(inputimg(:,:,1));
        G = double(inputimg(:,:,2));
        B = double(inputimg(:,:,3));
        undistimg(:,:,1) = round(interp2(R,u_d,v_d)); % R
        undistimg(:,:,2) = round(interp2(G,u_d,v_d)); % G
        undistimg(:,:,3) = round(interp2(B,u_d,v_d)); % B
    end
end 