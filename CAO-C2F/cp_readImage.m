
function [I1_gray, I2_gray, I1_rgb, I2_rgb, Filename1, Filename2, pathname] = cp_readImage(distort,inputpath,file1,file2)

%%
clc
close all
%% Read images needed to be registrated
disp('【0】Reading images ...')
filenameset = ['JPG' 'PNG' 'TIF' 'BMP' 'PEG' 'PPM' 'PGM'...
               'jpg' 'png' 'tif' 'bmp' 'peg' 'ppm' 'pgm' 'dat'];
ptemp = matlabpath;
ind = strfind(ptemp,';');
filepath = ptemp(ind(end)+1:end);

if nargin == 0
    distort = 0;
    [Filename1,pathname] = uigetfile([filepath '\*.*'],'pick an infrared image');
    if ~contains(filenameset,Filename1(end-2:end))
        error('Incorrect file format input!');
    end
    if sum( 'dat'== Filename1(end-2:end))==3
        [I1,I11,I2,pathname,Filename1,Filename2] = cp_readTempMat(pathname,Filename1);
    else
        I1 = imread([pathname Filename1]);
        [Filename2,pathname] = uigetfile([pathname '*' Filename1(end-3:end)],'pick a visible image');
        I2 = imread([pathname Filename2]);
    end
elseif nargin == 1
    [Filename1,pathname] = uigetfile([filepath '\*.*'],'pick an infrared image');
    if ~contains(filenameset,Filename1(end-2:end))
        error('Incorrect file format input!');
    end
    if sum( 'dat'== Filename1(end-2:end))==3
        [I1,I11,I2,pathname,Filename1,Filename2] = cp_readTempMat(pathname,Filename1);
    else
        I1 = imread([pathname Filename1]);
        [Filename2,pathname] = uigetfile([pathname '*' Filename1(end-3:end)],'pick a visible image');
        I2 = imread([pathname Filename2]);
    end
elseif nargin  == 4
    I1 = imread([inputpath file1]);
    I2 = imread([inputpath file2]);
    pathname = inputpath;
    Filename1 = file1;
    Filename2 = file2;
else
    error('No sufficient input parameters!');
end

path(path, pathname(1:end-1));

if distort
    if size(I1,1)~=768 && size(I1,1)>288 && size(I1,1)<1080 && size(I2,1)==1190
        I1(289:end,:,:)  = [];
        I2(1081:end,:,:) = [];
    elseif ( size(I1,1)==768 && size(I1,2)==576 )||( size(I1,1)==800 && size(I1,2)==600 )
        I1 = cp_undistortimg(I1);
        I2 = cp_undistortimg(I2);
        disp('   Undistort image successfully!');
    elseif size(I1,2)==768 && size(I1,1)==576
        I1 = cp_undistortimg(imrotate(I1,90));
        I2 = cp_undistortimg(I2);
        disp('   Undistort image successfully!');
    elseif size(I2,2)==768 && size(I2,1)==576
        I2 = cp_undistortimg(imrotate(I2,90));
        I1 = cp_undistortimg(I1);
        disp('   Undistort image successfully!');
    end
end
fprintf('\tsize of Image 1 is [%d×%d]\n    size of Image 2 is [%d×%d]\n',size(I1,1),size(I1,2),size(I2,1),size(I2,2));
if size(I1,3)==3 & size(I2,3)==3 
        I1_rgb = I1;
        I2_rgb = I2;
        I1_gray = double(rgb2gray(I1));
        I2_gray = double(rgb2gray(I2));
elseif size(I1,3)==3 & size(I2,3)==1 
        I1_gray = double(rgb2gray(I1));
        I2_gray = double(I2);
        I1_rgb = I1;
        I2_rgb = 0;
elseif size(I1,3)==1 & size(I2,3)==3 
        I1_gray = double(I1);
        I2_gray = double(rgb2gray(I2));
        I1_rgb = 0;
        I2_rgb = I2;
else
        I1_gray = double(I1);
        I2_gray = double(I2);
        I1_rgb = 0;
        I2_rgb = 0;
end
I1_gray = uint8((I1_gray-min(min(I1_gray)))./(max(max(I1_gray))-min(min(I1_gray))) * 255);
I2_gray = uint8((I2_gray-min(min(I2_gray)))./(max(max(I2_gray))-min(min(I2_gray))) * 255);

disp('【0】Completed to read image！')
