% Creates a dataset of support features to use for training a support
% classifier.
addpath('common/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

load(consts.splitsPath, 'trainNdxs');

% We need the labels just to check on where the floor is.
load(consts.datasetPath, 'namesToIds', 'names');
floorId = namesToIds('floor');

% Load the support labels.
load(consts.supportLabels, 'supportLabels');

%%
allFeatures = [];
allLabels = [];
allTrainIndics = [];

for ii = 1 : consts.numImages
  fprintf('Loading support features %d/%d.\n', ii, consts.numImages);
  if ~consts.useImages(ii)
    continue;
  end
  
  load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
  
  % Load ground-truth regions.
  imgRegionsTrue = get_regions(ii);
  
  % Load the region map we're actually using.
  imgRegionsActual = get_regions(ii, params);
  
  % Figure out which regions represent the floor. Since the segmentation
  % may be imperfect, only use those regions that are a majority floor.
  floorRegions = [];
  [regionIds, R] = get_region_ids(imgRegionsActual);
  floorMask = imgObjectLabels == floorId;
  for rr = 1  : R
    regionMask = imgRegionsActual == regionIds(rr);
    if nnz(regionMask & floorMask) / nnz(regionMask) > .5
      floorRegions = [floorRegions; regionIds(rr)];
    end
  end
  floorRegions(floorRegions == 0) = [];
  
  [features, supporteeIds, supporterIds] = load_support_features(ii, params);
  
  F = size(features, 1);
  
  % Select the ground truth.
  supportRelnsCur = supportLabels{ii};
  supportRelnsCur = get_ground_truth_support_relns(imgRegionsTrue, ...
    imgRegionsActual, supportRelnsCur, imgObjectLabels, floorId);
  
  % Now, figure out which of the features represents a support
  % relationship.
  supportVertRelnsForImage = supportRelnsCur(supportRelnsCur(:,3) == 1, 1:2);
  supportHorzRelnsForImage = supportRelnsCur(supportRelnsCur(:,3) == 2, 1:2);
  supportRelnsFromFeatures = [supporteeIds supporterIds];
  
  [~, ~, indVert] = intersect(supportVertRelnsForImage, supportRelnsFromFeatures, 'rows');
  [~, ~, indHorz] = intersect(supportHorzRelnsForImage, supportRelnsFromFeatures, 'rows');
  
  featureLabels = ones(F,1);
  featureLabels(indVert) = 2;
  featureLabels(indHorz) = 3;
  featureLabels(isin(supporteeIds, floorRegions)) = 4;
  
  allFeatures = [allFeatures; features];
  allLabels = [allLabels; featureLabels];
  allTrainIndics = [allTrainIndics; ones(F,1) * isin(ii, trainNdxs)];
end

allTrainIndics = logical(allTrainIndics);

%% Split into train and test.
trainData = allFeatures(allTrainIndics, :);
trainLabels = allLabels(allTrainIndics);

testData = allFeatures(~allTrainIndics, :);
testLabels = allLabels(~allTrainIndics);

%% Save to disk.
outFilename = sprintf(consts.supportDataset, params.regionSrc, ...
    params.seg.featureSet, params.stage);
save(outFilename, 'trainData', 'trainLabels', 'testData', 'testLabels', '-v7.3');
