% Sample code for detecting Harris corners, following 
% Brown et al, CVPR 2005
%        by Alyosha Efros
%
% function [x,y,strength]= harris(im,topn);
%
% Input: 
% im ... grayscale image
% topn ... how many features to detect
% 
% Ouput: x,y ... feature co-ordinates, strength==Harris corner strength
%
function [x,y,strength]= harris(im,topn);
%im = im2double(rgb2gray(imrgb));
g1 = fspecial('gaussian', 9,1);  % Gaussian with sigma_d 
g2 = fspecial('gaussian', 11,1.5); % Gaussian with sigma_i 
img1 = conv2(im,g1,'same');  % blur image with sigma_d
Ix = conv2(img1,[-1 0 1],'same');  % take x derivative 
Iy = conv2(img1,[-1;0;1],'same');  % take y derivative

% Compute elements of the Harris matrix H
%%% we can use blur instead of the summing window
Ix2 = conv2(Ix.*Ix,g2,'same');
Iy2 = conv2(Iy.*Iy,g2,'same');
IxIy = conv2(Ix.*Iy,g2,'same');
R = (Ix2.*Iy2 - IxIy.*IxIy) ... % det(H) 
    ./ (Ix2 + Iy2 + eps);       % trace(H) + epsilon

% don't want corners close to image border
R([1:20, end-21:end], :) = 0;
R(:,[1:20,end-21:end]) = 0;

% non-maxima supression within 3x3 windows
nonmax = inline('max(x)');
Rmax = colfilt(R,[3 3],'sliding',nonmax); % find neighbrhood max
Rnm = R.*(R == Rmax);  % supress non-max

% extract all interest points
[y,x,strength] = find(Rnm);

% return the top topn points only
[yy ii]  = sort(-strength);
nn       = min(topn,length(ii));
y        = y(ii(1:nn));
x        = x(ii(1:nn));
strength = strength(ii(1:nn));


% % show 'em
% figure,imagesc(im);
% colormap(gray);
% hold on;
% plot(x,y,'r.');
% hold off;
% 
