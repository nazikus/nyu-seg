% Determines the class and instance labels for each of the initial superpixel segments produced by
% the watershed segmentation.

addpath('common/');
addpath('segmentation/');
Consts;
Params;

OVERWRITE = params.overwrite_labels;

%%
fprintf('\nRunning regions2labels on superpixels from Watershed:\n');
global ii_;

for ii_ = 1 : consts.numImages
  if ~consts.useImages(ii_)
    continue;
  end
  
  fprintf('running regions2labels %d/%d... ', ii_, consts.numImages);
  outFilename = sprintf(consts.objectLabelsSegFilename, ii_);
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('skipping (exists), overwrite=false.\n');
    continue;
  end
  
  fprintf('Running regions2labels (%d/%d)\r', ii_, consts.numImages);
  
  
  load(sprintf(consts.objectLabelsFilename, ii_), 'imgObjectLabels');
  load(sprintf(consts.instanceLabelsFilename, ii_), 'imgInstanceLabels');

  load(sprintf(consts.watershedFilename, ii_), 'boundaryInfo');
  
  [instanceMasks, instanceLabels] = get_instance_masks(imgObjectLabels, imgInstanceLabels);
  %MARKER region2labels
  [classLabels, instanceLabels, intersectionPcnt] = ...
      regions2labels(boundaryInfo.imgRegions, instanceMasks, instanceLabels);
    
  save(outFilename, 'classLabels', 'instanceLabels', 'intersectionPcnt');
  fprintf(' Done.\n');
end

clearvars -global ii_;

fprintf('\n');
fprintf('===============================\n');
fprintf('Finished running region2labels \n');
fprintf('===============================\n\n');