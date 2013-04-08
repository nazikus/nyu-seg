% Create the regions from the ground truth (manually annotated) labels.
addpath('common/');
Consts;

if ~exist(consts.imageRegionsDir, 'dir')
  mkdir(consts.imageRegionsDir);
end

%%
for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  fprintf('Extracting regions from labels %d/%d.\n', ii, consts.numImages);
  
  load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
  load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
  imgRegions = get_regions_from_labels(imgObjectLabels, imgInstanceLabels);
  
  % Chuck any regions that are not big enough.
  [~, R] = get_region_ids(imgRegions);
  
  for rr = 1 : R
    regionMask = imgRegions == rr;
    if nnz(regionMask) < consts.MIN_PIXELS_PER_REGION
      fprintf('Discarding a region with less than %d pixels.\n', consts.MIN_PIXELS_PER_REGION);
      imgRegions(regionMask) = 0;
    end
  end
  
  % Now, go back and remap the labels to avoid situations in which we have
  % regions [1 2 4].
  [regionIds, R] = get_region_ids(imgRegions);
  
  for rr = 1 : R
    regionMask = imgRegions == regionIds(rr);
    imgRegions(regionMask) = rr;
  end
  
  [regionIds, R] = get_region_ids(imgRegions);
  assert(max(imgRegions(:)) == R);
  assert(max(imgRegions(:)) == numel(regionIds));
  
  save(sprintf(consts.imageRegionsFilename, ii), 'imgRegions');
end
