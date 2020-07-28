 
function [Iaffine,affmat] = cp_getAffine(I1,I2,p1,p2)
if size(p1,1) == 3 
    affmat = fitgeotrans(p1,p2,'affine');
    disp('     Affine transformation applied!');
    Iaffine = imwarp(I1,affmat,'OutputView',imref2d(size(I2)));
elseif size(p1,1) >= 4
    affmat = fitgeotrans(p1,p2,'projective');
    disp('     Projective transformation applied!');
    Iaffine = imwarp(I1,affmat,'OutputView',imref2d(size(I2)));
elseif size(p1,1) == 2 
    affmat = fitgeotrans(p1,p2,'nonreflectivesimilarity');
    disp('     Nonreflective similarity transformation applied!');
    Iaffine = imwarp(I1,affmat,'OutputView',imref2d(size(I2)));
else
    error('Transformation Failed! No sufficient Matches!!');
end
