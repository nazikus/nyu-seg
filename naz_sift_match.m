clear; clc; close all;

addpath(genpath('.\'));
ind_ = @(A,r,c) A(r,c); 

Consts; Params;
params.debug_visible = 'on';   % doesn't work because seg2framents.m loads Params.m again

cosegImgNdx = {[1 2]; [3 4]; [5 6 7]; [8]; [9]; [10 11 12] };

for i = cosegImgNdx{6}
    load(sprintf(consts.imageRgbFilename,  i), 'imgRgb');
    load(sprintf(consts.watershedFilename, i), 'boundaryInfo');
    %MARKER round!
    juncs = round(boundaryInfo.junctions.position);  % round means top-right pixel from junction point
    imgRgb = double(imgRgb)/255; % im2double
    
    switch i
        case cosegImgNdx{6}(1)
            imargb = imgRgb;
            xa = juncs(:,1);
            ya = juncs(:,2);
            
        case cosegImgNdx{6}(2)
            imbrgb = imgRgb;
            xb = juncs(:,1);
            yb = juncs(:,2);
            
        case cosegImgNdx{6}(3)
            imcrgb = imgRgb;
            xc = juncs(:,1);
            yc = juncs(:,2);
    end    
end
clear imgRgb boundaryInfo;

%MARKER just to avoid error at dist2 function (line 76)
% ====================================================
% common_len = min([length(xa),length(xb),length(xc)]);
% %aNxs = ind_(randperm(length(xa)),1,common_len);
% Ndx = randperm(length(xa));
% Ndx = Ndx(1:common_len);
% xa = xa(Ndx);
% ya = ya(Ndx);
% 
% Ndx = randperm(length(xb));
% Ndx = Ndx(1:common_len);
% xb = xb(Ndx);
% yb = yb(Ndx);
% 
% Ndx = randperm(length(xc));
% Ndx = Ndx(1:common_len);
% xc = xc(Ndx);
% yc = yc(Ndx);
% ====================================================
%END

ima = rgb2gray(imargb);
imb = rgb2gray(imbrgb);
imc = rgb2gray(imcrgb);

% show images
figure(100); clf;
subplot(1,3,1); imagesc(imargb); axis image; axis off; title('Image a');
subplot(1,3,2); imagesc(imbrgb); axis image; axis off; title('Image b');
subplot(1,3,3); imagesc(imcrgb); axis image; axis off; title('Image c');

% show detected points
figure(1); clf;
imagesc(imargb); axis image; hold on;
plot(xa,ya,'+y');

% show all points
figure(2); clf;
imagesc(imbrgb); axis image; hold on;
plot(xb,yb,'+y');

% show all points
figure(3); clf;
imagesc(imcrgb); axis image; hold on;
plot(xc,yc,'+y');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract descriptors (heavily blurred 21xb1 patches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Da,xa,ya] = ext_desc(ima,xa,ya);
[Db,xb,yb] = ext_desc(imb,xb,yb);
[Dc,xc,yc] = ext_desc(imc,xc,yc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute tentative matches between image 1 (a) and 2 (b) 
% by matching local features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rt        = 0.4;              % 1NN/2NN distance ratio threshold (between 0 and 1)
D2        = dist2(Da',Db'); % compute pair-wise distances between descriptors
[Y,I]     = sort(D2,2);     % sort distances
rr        = Y(:,1)./Y(:,2); % compute D. Lowes' 1nn/2nn ratio test
inD12     = find(rr<rt);   % take only points with a 1nn/2nn ratio below 0.8
I         = I(inD12);       % select matched points
xat       = xa(inD12);
yat       = ya(inD12);
xbt       = xb(I);
ybt       = yb(I);

% show all tentative matches
figure(1); clf;
imagesc(imargb); hold on;
plot(xat,yat,'+g');
hl = line([xat; xbt],[yat; ybt],'color','y');
title('Tentative correspondences');
axis off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robustly fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the inlier threshold (in noramlized image co-ordinates)
t              = 0.3; 
[Hab, inliers] = ransacfithomography([xat; yat], [xbt; ybt], t);

% show inliers
figure(4); clf;
imagesc(imargb); hold on;
hl = line([xat(inliers); xbt(inliers)],[yat(inliers); ybt(inliers)]);
set(hl,'color','y');
plot(xat(inliers),yat(inliers),'+g');
title('Inliers');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vgg_gui_H(imc,imb,Hcb);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warp and composite images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf;
bbox=[-400 1200 -200 700]   % image space for mosaic
% warp image b to mosaic image using an identity homogrpahy
% Image b is chosen as the reference frame
iwb = vgg_warp_H(imbrgb, eye(3), 'linear', bbox);
imshow(iwb); axis image;

% warp image 1 to the reference mosaic frame (image 2) 
figure(6); clf;
iwa = vgg_warp_H(imargb, Hab, 'linear', bbox);  % warp image a to the mosaic image
imshow(iwa); axis image;
imagesc(double(max(iwb,iwa))); % combine images into a common mosaic (take maximum value of the two images)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate homography between images 3 and 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on the code above, write code to:
% 1. Compute tentative matches between images 2 and 3 
% 2. Robustly fit homography Hcb
% 3. Re-estimate homography from inliers


% ---- 

rt        = 0.5;              % 1NN/2NN distance ratio threshold (between 0 and 1)
D2        = dist2(Dc',Db'); % compute pair-wise distances between descriptors
[Y,I]     = sort(D2,2);     % sort distances
rr        = Y(:,1)./Y(:,2); % compute D. Lowes' 1nn/2nn ratio test
inD12     = find(rr<rt);   % take only points with a 1nn/2nn ratio below 0.8
I         = I(inD12);       % select matched points
xct       = xc(inD12);
yct       = yc(inD12);
xbt       = xb(I);
ybt       = yb(I);

% show all tentative matches
figure(1); clf;
imagesc(imargb); hold on;
plot(xat,yat,'+g');
hl = line([xct; xbt],[yct; ybt],'color','y');
title('Tentative correspondences');
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robustly fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the inlier threshold (in noramlized image co-ordinates)
t              = 0.3; 
[Hcb, inliers] = ransacfithomography([xct; yct], [xbt; ybt], t);

% show inliers
figure(4); clf;
imagesc(imargb); hold on;
hl = line([xct(inliers); xbt(inliers)],[yct(inliers); ybt(inliers)]);
set(hl,'color','y');
plot(xct(inliers),yct(inliers),'+g');
title('Inliers');

figure(5); clf;
bbox=[-400 1200 -200 700]   % image space for mosaic
% warp image b to mosaic image using an identity homogrpahy
% Image b is chosen as the reference frame
iwb = vgg_warp_H(imbrgb, eye(3), 'linear', bbox);
imshow(iwb); axis image;

% warp image 1 to the reference mosaic frame (image 2) 
figure(6); clf;
iwc = vgg_warp_H(imcrgb, Hcb, 'linear', bbox);  % warp image a to the mosaic image
imshow(iwc); axis image;
imagesc(double(max(iwb,iwc))); % combine images into a common mosaic (take maximum value of the two images)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final warping and compositing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7); clf;
iwc = vgg_warp_H(imcrgb, Hcb, 'linear', bbox);  % warp image c to mosaic image
imshow(iwc); axis image;

figure(8); clf;
imagesc(max(iwc,double(max(iwb,iwa)))); % combine images into a common mosaic
axis image; axis off;



