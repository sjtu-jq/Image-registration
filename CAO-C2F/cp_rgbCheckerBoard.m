function cb = cp_rgbCheckerBoard(I1aff,I2rgb,unit)
%CP_GRAYFUSION obtain the rgb fusion image with checkerboard format
cb = I2rgb;
I1gray = rgb2gray(I1aff);
[v,u] = find(I1gray > 0);
v = sort(v); u = sort(u);
umin = u(1); umax = u(end);
vmin = v(1); vmax = v(end);
rowgap = ceil( (vmax-vmin)/unit);
colgap = ceil( (umax-umin)/unit);
for i = 1:unit
    for j = (2-mod(i,2)):2:unit
        cb(vmin+(i-1)*rowgap:min(vmax,vmin+i*rowgap),umin+(j-1)*colgap:min(umax,umin+j*colgap),:) =...
     I1aff(vmin+(i-1)*rowgap:min(vmax,vmin+i*rowgap),umin+(j-1)*colgap:min(umax,umin+j*colgap),:);
    end
end
[v,u] = find(cb(:,:,1) == 0);
for i = 1:length(v)
    if cb(v(i),u(i),2) ~= 0 || cb(v(i),u(i),3) ~= 0 
        continue;
    end
    cb(v(i),u(i),:) = I2rgb(v(i),u(i),:);
end

