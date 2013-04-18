% function [D,xD,yD] = ext_desc(Im,x,y)
%
% Extract patch descriptors around detected local features 
% Descriptors are 21x21 grayscale pixel patches blurred 
% by a (fairly large) Gaussian with sigma=3 pixels.
%
% Based on 
% “Multi-Image Matching using Multi-Scale Oriented Patches” by Brown et
% al., CVPR'2005
%
% Input: 
% Im  ... grayscale image
% x,y ... (npoints x 1) image co-ordinates of detected features
%
% Output: 
% D... (441 x npoints) matrix where each column is one vectorized
% descriptor of dimension 441 (21*21 pixels).
% xD,yD ... co-ordinates of each extracted patch.
%
% Josef.Sivic@ens.fr
%
function [D,xD,yD] = ext_desc(Im,x,y)
%

% blur images
sigma = 3;
B     = fspecial('gaussian',6*sigma+1,sigma);
Imb   = imfilter(Im,B,'same');


%figure(1); clf; imagesc(Imb); colormap gray;
w     = 10;
j = 1;
for i = 1:length(x)
    xi = round(x(i));
    yi = round(y(i));
    minx = xi-w;
    maxx = xi+w;
    miny = yi-w;
    maxy = yi+w;
    
    % patches reaching out of image bounds are ignored
    if min(minx,miny)<=0 || maxx>size(Imb,2) || maxy>size(Imb,1)
        continue;
    end;
    
    W = Imb(miny:maxy,minx:maxx);
    D(:,j) = W(:);
    xD(j)  = xi;
    yD(j)  = yi;
    j=j+1;    
end;
