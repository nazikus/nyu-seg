% Visualizes all of the images.
addpath('common/');
Consts;
Params;


for ii = 5 %1 : consts.numImages
  
  % Load the RGB image.
  load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
  
  % Load the Depth image.
  load(sprintf(consts.imageDepthFilename, ii), 'imgDepth');
  load(sprintf(consts.imageDepthRawFilename, ii), 'imgDepthRaw');

  % Load the Labels.
  imgRegions = get_regions(ii);
  imgRegionsS = get_regions(ii,params);
  
  
  sfigure(1);
%   subplot(1,3,1);
  imshow(imgRgb,[]);
  title('RGB');
  axis off;
  axis equal;
  
%   subplot(1,3,2);
%   imagesc(imgDepth);
%   title('Depth');
%   axis off;
%   axis equal;
  
%   subplot(1,3,2);
%   vals = randperm(max(imgRegions(:)));
%   vis_regions(imgRegions, vals);
%   title('Instances');
%   axis off;
%   axis equal;
% 
%   subplot(1,3,3);
%   vals = randperm(max(imgRegionsS(:)));
%   vis_regions(imgRegionsS, vals);
%   title('Segmented');
%   axis off;
%   axis equal;
  
  set(gcf, 'Name', sprintf('Image %d', ii));
  %pause;
end