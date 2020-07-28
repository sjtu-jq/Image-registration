function [cout, angle, edge_map]=cp_cornerDetection(varargin)

[I,C,T_angle,sig,H,L,Endpoint,Gap_size,maxlength,rflag] = parse_inputs(varargin{:});

BW=edge(I,'canny',[L,H]);  % Detect edges

[curve,curve_start,curve_end,curve_mode,curve_num,edge_map]=extract_curve(BW,Gap_size);  % Extract curves

[cout,angle]=get_corner(curve,curve_start,curve_end,curve_mode,curve_num,BW,sig,Endpoint,C,T_angle,maxlength,rflag); % Detect corners

function [curve,curve_start,curve_end,curve_mode,cur_num,edge_map]=extract_curve(BW,Gap_size)

%   Function to extract curves from binary edge map, if the endpoint of a
%   contour is nearly connected to another endpoint, fill the gap and continue
%   the extraction. The default gap size is 1 pixles.

[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1);
cur_num=0;
% 提取边缘图中所有的连续轮廓并存为元胞数组：
while size(r,1)>0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [~,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1; % find continous point
        cur=[cur;point];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [~,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
        
    if size(cur,1)>(size(BW,1)+size(BW,2))/60
        cur_num=cur_num+1;
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end
% 将轮廓元胞数组中的起始点存起来并判定是否为闭合轮廓
for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
        (curve_start(i,2)-curve_end(i,2))^2 <= 4 %32
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1; % 还原边缘图像
end
edge_map = BW_edge;


function [cout,angle]=get_corner(curve,curve_start,curve_end,curve_mode,curve_num,BW,sig,Endpoint,C,T_angle,maxlength,rflag)

corner_num=0;
cout=[];
angle=[];

GaussianDieOff = .0001; 
pw = 1:30; 
ssq = sig*sig;
width = max( find( exp(-(pw.*pw)/(2*ssq)) > GaussianDieOff) );
if isempty(width)
    width = 1;  
end
t = (-width:width); %高斯窗 width=12
gau = exp(-(t.*t)/(2*ssq)); 
gau=gau/sum(gau);
sigmaLs =  maxlength;
warning off
% beforegaus = zeros(size(BW));
% aftergaus = zeros(size(BW));
parfor i=1:curve_num
    x=curve{i}(:,1);
    y=curve{i}(:,2);
%     for k=1:length(x)
%        beforegaus(x(k),y(k)) = 1; 
%     end
    W=width;
    L=length(x);
    if L>W
        % Calculate curvature
        if curve_mode(i,:)=='loop'
            xL=[x(L-W+1:L);x;x(1:W)];
            yL=[y(L-W+1:L);y;y(1:W)];
        else
            xL=[ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];% 总长度2w+L
            yL=[ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end
       
        xx=conv(xL,gau);
        xx=xx(W+1:L+3*W);
        yy=conv(yL,gau);
        yy=yy(W+1:L+3*W);
%         aftergaus(:)=1; % show curves
%         for k=width:length(x)-width+1
%             if xx(k) < 1 || yy(k) < 1
%                 continue;
%             end
%            aftergaus(round(xx(k)),round(yy(k))) = 0; 
%         end
%         figure,imshow(~aftergaus);
        Xu=[xx(2)-xx(1) ; (xx(3:L+2*W)-xx(1:L+2*W-2))/2 ; xx(L+2*W)-xx(L+2*W-1)];
        Yu=[yy(2)-yy(1) ; (yy(3:L+2*W)-yy(1:L+2*W-2))/2 ; yy(L+2*W)-yy(L+2*W-1)];
        Xuu=[Xu(2)-Xu(1) ; (Xu(3:L+2*W)-Xu(1:L+2*W-2))/2 ; Xu(L+2*W)-Xu(L+2*W-1)];
        Yuu=[Yu(2)-Yu(1) ; (Yu(3:L+2*W)-Yu(1:L+2*W-2))/2 ; Yu(L+2*W)-Yu(L+2*W-1)];
        K=abs((Xu.*Yuu-Xuu.*Yu)./((Xu.*Xu+Yu.*Yu).^1.5));
        K=ceil(K*100)/100;
               
        % Find curvature local maxima as corner candidates
        extremum=[];
        N=size(K,1);
        n=0;
        Search=1;
        
        for j=1:N-1
            if (K(j+1)-K(j))*Search>0
                n=n+1;
                extremum(n)=j;  % In extremum, odd points is minima and even points is maxima
                Search=-Search;
            end
        end
        if mod(size(extremum,2),2)==0
            n=n+1;
            extremum(n)=N;
        end
    
        n=size(extremum,2);
        flag=ones(size(extremum));
        lambda_LR=[];
        % Compare with adaptive local threshold to remove round corners
        for k=2:2:n
            [~,index1]=min(K(extremum(k):-1:extremum(k-1)));
            [~,index2]=min(K(extremum(k):extremum(k+1)));
            ROS=K(extremum(k)-index1+1:extremum(k)+index2-1);
            K_thre = C*mean(ROS);
            if K(extremum(k))<K_thre
                flag(k)=0;
            else
                lambda_LR=[lambda_LR; extremum(k) index1 index2];
            end
        end
        extremum=extremum(2:2:n);
        flag=flag(2:2:n);
        extremum=extremum(find(flag==1));
        
        % Check corner angle to remove false corners due to boundary noise and trivial details
        smoothed_curve=[xx,yy];
        n=size(extremum,2);
        flag=ones(size(extremum)); 
        for j=1:n
            if j==1 & j==n
                ang=curve_tangent(smoothed_curve(1:L+2*W,:),extremum(j));
            elseif j==1 
                ang=curve_tangent(smoothed_curve(1:extremum(j+1),:),extremum(j));
            elseif j==n
                ang=curve_tangent(smoothed_curve(extremum(j-1):L+2*W,:),extremum(j)-extremum(j-1)+1);
            else
                ang=curve_tangent(smoothed_curve(extremum(j-1):extremum(j+1),:),extremum(j)-extremum(j-1)+1);
            end     
            if ang>T_angle & ang<(360-T_angle)
                flag(j)=0;  
            end
        end
        extremum=extremum(find(flag~=0));
        lambda_LR=lambda_LR(find(flag~=0),:);
            
        extremum=extremum-W;
        true_corner = find(extremum>0 & extremum<=L);
        extremum=extremum(true_corner);
        lambda_LR = lambda_LR(true_corner,:);
        n=size(extremum,2);     
        for j=1:n     
            cout = [cout; curve{i}(extremum(j),:)];
            if rflag                
%**********    CAO: adaptive  parameters algorithm  *********% 
                xcor = curve{i}(extremum(j),2);
                ycor = curve{i}(extremum(j),1);
                retail = size(curve{i},1);
                lengthL = min([lambda_LR(j,2) extremum(j)]);            
                lengthR = min([lambda_LR(j,3) retail-extremum(j)+1]);
                coefL = exp(-((0:1:lengthL-1) / lengthL).^2/2);
                coefL = coefL / sum(coefL);
                coefR = exp(-((lengthR-1:-1:0) / lengthR).^2/2);
                coefR = coefR / sum(coefR);
                xL= sum(coefL' .* curve{i}(extremum(j)+1-lengthL:extremum(j),2));
                yL= sum(coefL' .* curve{i}(extremum(j)+1-lengthL:extremum(j),1));
                xR = sum(coefR' .* curve{i}(extremum(j):extremum(j)+lengthR-1,2));
                yR = sum(coefR' .* curve{i}(extremum(j):extremum(j)+lengthR-1,1));
                vL = [xL-xcor yL-ycor]; vR = [xR-xcor yR-ycor];
                vm = vL / norm(vL) + vR / norm(vR);
%*************   Assign Main Orientation  **********%
                deltax = vm(1);
                deltay = vm(2);
                if isnan(atan(deltay/deltax)) % no solution
                    orientation1 = 0;
                elseif deltay >= 0 & deltax >= 0 % quadrant I
                    orientation1 = atan(deltay / deltax);
                elseif deltay < 0 & deltax >= 0  % quadrant IV
                    orientation1 = atan(deltay / deltax) + 2*pi;
                elseif (deltay >= 0 & deltax < 0) | (deltay < 0 & deltax < 0 ) % quadrant II & III
                    orientation1 = atan(deltay / deltax) + pi;
                end
                angle = [angle; orientation1];
            end
        end
    end
end
% figure,subplot(121),imshow(~beforegaus);subplot(122),imshow(~aftergaus);

% Add Endpoints
if Endpoint
    for i=1:curve_num
        retail = size(curve{i},1);
        if retail>0 & curve_mode(i,:)=='line'
            % Start point compare with detected corners
            compare_corner=cout-ones(size(cout,1),1)*curve_start(i,:);
            compare_corner=compare_corner.^2;
            compare_corner=compare_corner(:,1)+compare_corner(:,2);
            if min(compare_corner)>25       % Add end points far from detected corners 
                cout = [cout; curve_start(i,:)];
                if rflag
                    xx = curve{i}(1,2);
                    yy = curve{i}(1,1);
                    coef2 = exp(-((min(retail,maxlength)-1:-1:0) / sigmaLs).^2/2);
                    coef2 = coef2 / sum(coef2);
                    xL = sum(coef2' .*curve{i}(1:min(retail,maxlength),2));
                    yL = sum(coef2' .*curve{i}(1:min(retail,maxlength),1));
                    deltax = xL - xx;
                    deltay = yL - yy;
                    if isnan(atan(deltay/deltax)) % no solution
                        orientation2 = 0;
                    elseif deltay >= 0 & deltax >= 0 % quadrant I
                        orientation2 = atan(deltay / deltax);
                    elseif deltay < 0 & deltax >= 0  % quadrant IV
                        orientation2 = atan(deltay / deltax) + 2*pi;
                    elseif (deltay >= 0 & deltax < 0) | (deltay < 0 & deltax < 0 ) % quadrant II & III
                        orientation2 = atan(deltay / deltax) + pi;
                    end
                    angle = [angle;orientation2];
                end
            end
            
            % End point compare with detected corners
            compare_corner=cout-ones(size(cout,1),1)*curve_end(i,:);
            compare_corner=compare_corner.^2;
            compare_corner=compare_corner(:,1)+compare_corner(:,2);
            if min(compare_corner)>25
                cout = [cout; curve_end(i,:)]; 

                if rflag
                    xx = curve{i}(end,2);
                    yy = curve{i}(end,1);
                    coef1 = exp(-((0:1:retail-max(1,retail-maxlength+1)) / sigmaLs).^2/2);
                    coef1 = coef1 / sum(coef1);
                    xL = sum(coef1' .*curve{i}(max(1,retail-maxlength+1):retail,2));
                    yL = sum(coef1' .*curve{i}(max(1,retail-maxlength+1):retail,1));
                    deltax = xL - xx;
                    deltay = yL - yy;
                    if isnan(atan(deltay/deltax)) % no solution
                        orientation3 = 0;
                    elseif deltay >= 0 & deltax >= 0 % quadrant I
                        orientation3 = atan(deltay / deltax);
                    elseif deltay < 0 & deltax >= 0  % quadrant IV
                        orientation3 = atan(deltay / deltax) + 2*pi;
                    elseif (deltay >= 0 & deltax < 0) | (deltay < 0 & deltax < 0 ) % quadrant II & III
                        orientation3 = atan(deltay / deltax) + pi;
                    end
                    angle = [angle;orientation3];
                end
            end
        end
    end
end

function ang=curve_tangent(cur,center)
for i=1:2
    if i==1
        curve=cur(center:-1:1,:);
    else
        curve=cur(center:size(cur,1),:);
    end
    L=size(curve,1);
    
    if L>3
        if sum(curve(1,:)~=curve(L,:))~=0
            M=ceil(L/2);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M,1);
            y2=curve(M,2);
            x3=curve(L,1);
            y3=curve(L,2);
        else
            M1=ceil(L/3);
            M2=ceil(2*L/3);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M1,1);
            y2=curve(M1,2);
            x3=curve(M2,1);
            y3=curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2*(-y1*x2^2+y3*x2^2-y3*y1^2-y3*x1^2-y2*y3^2+x3^2*y1+y2*y1^2-y2*x3^2-y2^2*y1+y2*x1^2+y3^2*y1+y2^2*y3)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            y0 = -1/2*(x1^2*x2-x1^2*x3+y1^2*x2-y1^2*x3+x1*x3^2-x1*x2^2-x3^2*x2-y3^2*x2+x3*y2^2+x1*y3^2-x1*y2^2+x3*x2^2)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_direction=angle(complex(x0-x1,y0-y1));
            adjacent_direction=angle(complex(x2-x1,y2-y1));
            tangent_direction=sign(sin(adjacent_direction-radius_direction))*pi/2+radius_direction;
        end
    
    else % very short line
        tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    direction(i)=tangent_direction*180/pi;
end
ang=abs(direction(1)-direction(2));


% 用矩形框标记角点

function img1=mark(img,x,y,w)

[M,N,C]=size(img);
img1=img;

if isa(img,'logical')
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
else
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<128)*255;
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
end

% 分解输入的实际参数
function [I,C,T_angle,sig,H,L,Endpoint,Gap_size,maxlength,rflag] = parse_inputs(varargin);

narginchk(0,10); % 判定输入参数数目是否在0~10之间
Para=[1.5,170,3,0.2,0,1,1,2,0]; %Default experience value;
if nargin>=2
    I=varargin{1};
    for i=2:nargin
        if size(varargin{i},1)>0
            Para(i-1)=varargin{i};
        end
    end
end

if nargin==1
    I=varargin{1};
end
    
if nargin==0 | size(I,1)==0
    [fname,dire]=uigetfile('*.bmp;*.jpg;*.gif','Open the image to be detected');
    I=imread([dire,fname]);
end

C=Para(1);
T_angle=Para(2);
sig=Para(3);
H=Para(4);
L=Para(5);
Endpoint=Para(6);
Gap_size=Para(7);
maxlength = Para(8);
rflag = Para(9);

%   CORNER Find corners in intensity image. 
%   
%       CORNER works by the following step:
%       1.	Apply the Canny edge detector to the gray level image and obtain a
%       binary edge-map.
%       2.	Extract the edge contours from the edge-map, fill the gaps in the
%       contours.
%       3.	Compute curvature at a low scale for each contour to retain all
%       true corners.
%       4.	All of the curvature local maxima are considered as corner
%       candidates, then rounded corners and false corners due to boundary
%       noise and details were eliminated.
%       5.  End points of line mode curve were added as corner, if they are not
%       close to the above detected corners.
%
%       Syntax :    
%       [cout,marked_img]=corner(I,C,T_angle,sig,H,L,Endpiont,Gap_size)
%
%       Input :
%       I -  the input image, it could be gray, color or binary image. If I is
%           empty([]), input image can be get from a open file dialog box.
%       C -  denotes the minimum ratio of major axis to minor axis of an ellipse, 
%           whose vertex could be detected as a corner by proposed detector.  
%           The default value is 1.5.
%       T_angle -  denotes the maximum obtuse angle that a corner can have when 
%           it is detected as a true corner, default value is 162.
%       Sig -  denotes the standard deviation of the Gaussian filter when
%           computeing curvature. The default sig is 3.
%       H,L -  high and low threshold of Canny edge detector. The default value
%           is 0.35 and 0.
%       Endpoint -  a flag to control whether add the end points of a curve
%           as corner, 1 means Yes and 0 means No. The default value is 1.
%       Gap_size -  a paremeter use to fill the gaps in the contours, the gap
%           not more than gap_size were filled in this stage. The default 
%           Gap_size is 1 pixels.
%
%       Output :
%       cout -  a position pair list of detected corners in the input image.
%       marked_image -  image with detected corner marked.
%
%       Examples
%       -------
%       I = imread('alumgrns.tif');
%       cout = corner(I,[],[],[],0.2);
%
%       [cout, marked_image] = corner;
%
%       cout = corner([],1.6,155);
%
%
%   Composed by He Xiaochen 
%   HKU EEE Dept. ITSR, Apr. 2005
%
%   Algorithm is derived from :
%       X.C. He and N.H.C. Yung, Curvature Scale Space Corner Detector with  
%       Adaptive Threshold and Dynamic Region of Support, Proceedings of the
%       17th International Conference on Pattern Recognition, 2:791-794, August 2004.
%   Improved algorithm is included in :
%   	X.C. He and N.H.C. Yung, corner detector based on global and local curvature properties, 
%       Optical Engineering, 47(5), pp: 057008, 2008.
%  