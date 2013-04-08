% Extracts and saves all of the region-to-structure class features using the ground truth (manually annotated) regions.
addpath('common/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

addpath(consts.spamsPath);


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
  fprintf('Extracting region-to-structure-class features %d/%d.\n', ii, consts.numImages);
  
  outFilename = sprintf(consts.structureFeaturesFilename, ...
    params.regionSrc, params.seg.featureSet, params.stage, ii);
  
  if exist(outFilename, 'file') && ~OVERWRITE
    continue;
  end

  regionFeatures = extract_region_to_structure_classes_features(ii, params);
  assert(~any(isnan(regionFeatures(:))));
  assert(~any(isinf(regionFeatures(:))));
  save(outFilename, 'regionFeatures');
end
