% Performs an initial watershed segmentation on each RGBD image.
addpath(genpath('iccv07Final'));
addpath('segmentation/');

Consts;
Params;

OVERWRITE = params.overwrite_watershed;

%
if ~exist(consts.watershedDir, 'dir')
  mkdir(consts.watershedDir);
end
global ii_;

for ii_ = 1 : consts.numImages
  
  if ~consts.useImages(ii_)
    continue;
  end
  
  fprintf('running watershed %d/%d... ', ii_, consts.numImages);
  outFilename = sprintf(consts.watershedFilename, ii_);
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('skipping (exists), OVERWRITE=false.\n');
   continue;
  end

  load(sprintf(consts.imageRgbFilename, ii_), 'imgRgb');  
  load(sprintf(consts.planeDataFilename, ii_), 'planeData');
  
  [boundaryInfo, pbAll] = im2superpixels(imgRgb, double(planeData.planeMap));
  save(outFilename, 'boundaryInfo', 'pbAll');
  fprintf('done.\n');

end
clearvars -global ii_;

fprintf(' ================================================= \n');
fprintf('Finished initial watershed segmentation.\n');
fprintf(' ================================================= \n\n');
