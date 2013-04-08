addpath('common/');
% addpath('segmentation/'); % Remove this!
Consts;

OVERWRITE = 0;

%%
for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  
  fprintf('Labeling watershed segments (%d/%d).\n', ii, consts.numImages);
  outFilename = sprintf(consts.objectLabelsSegFilename, ii);
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('Skipping file %d/%d since it already exists.\n', ii, consts.numImages);
    continue;
  end

  load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
  load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');

  load(sprintf(consts.watershedFilename, ii), 'boundaryInfo');

  [classLabels, instanceLabels, intersectionPcnt] = ...
      get_labels_from_regions(boundaryInfo.imgRegions, imgObjectLabels, ...
        imgInstanceLabels);
  
  save(outFilename, 'classLabels', 'instanceLabels', 'intersectionPcnt');
end

fprintf('Finished labeling watershed segmentation.\n');
