% Extracts and saves all of the region-to-structure class features.
addpath('common/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;

addpath(consts.spamsPath);

% Whether or not to overwrite the structure class features file if it
% already exists on disk.
OVERWRITE = params.overwrite_struct;

if ~exist(consts.structureFeaturesDir, 'dir')
  mkdir(consts.structureFeaturesDir);
end

%%
RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));

for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  fprintf('Extracting region-to-structure-class features %d/%d... ', ii, consts.numImages);
  
  outFilename = sprintf(consts.structureFeaturesFilename, ...
    params.regionSrc, params.seg.featureSet, params.stage, ii);
  
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('skipping (exists), overwrite=false.\n');
    continue;
  end

  regionFeatures = extract_region_to_structure_classes_features(ii, params);
  assert(~any(isnan(regionFeatures(:))));
  assert(~any(isinf(regionFeatures(:))));
  save(outFilename, 'regionFeatures');
  fprintf('Done.\n');
end
