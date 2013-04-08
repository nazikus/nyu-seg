% Trains several stages of boundary classifiers
addpath('common/');
addpath('segmentation/');
addpath(genpath('iccv07Final'));
Consts;
Params;

params.seg.featureSet = consts.BFT_RGBD;

% Load the train/test split.
load(consts.splitsPath, 'trainNdxs');

if ~exist(consts.boundaryFeaturesDir, 'dir')
  mkdir(consts.boundaryFeaturesDir);
end

for stage = 1 : params.seg.numStages  
  %MARKER extract boundary classifier features and labels  
  extract_boundary_classifier_features_and_labels(stage, params);

  %MARKER Create the boundary-classification dataset.
  % datasetFilename = sprintf(consts.boundaryFeaturesDataset, params.seg.featureSet, stage);
  datasetFilename = sprintf(consts.boundaryFeaturesDataset, params.seg.featureSet, stage);
  datasetFilename = sprintf('%s_s%d.mat', datasetFilename(1:end-4), length(consts.useNdx));
  
  if ~exist(datasetFilename, 'file') || params.overwrite_train
    [trainData, testData, trainLabels, testLabels] = ...
        create_boundary_classifier_dataset(stage, trainNdxs, params.seg.featureSet);
    fprintf('Saving dataset...');
    save(datasetFilename, 'trainData', 'trainLabels', 'testData', 'testLabels', '-v7.3');
    fprintf('DONE\n');
  else
    fprintf('Loading the boundary-classification dataset.\n');
    load(datasetFilename, 'trainData', 'trainLabels', 'testData', 'testLabels');
  end

  %MARKER Train the boundary classifier.
  % boundaryClassifierFilename = sprintf(consts.boundaryClassifierFilename, params.seg.featureSet, stage);
  boundaryClassifierFilename = sprintf(consts.boundaryClassifierFilename, params.seg.featureSet, stage);
  boundaryClassifierFilename = sprintf('%s_s%d.mat', boundaryClassifierFilename(1:end-4), length(consts.useNdx) );
  
  if ~exist(boundaryClassifierFilename, 'file') || params.overwrite_train
    classifier = train_boundary_classifier_dt(stage, trainData, trainLabels, testData, testLabels, params);
    save(boundaryClassifierFilename, 'classifier');
  else
    fprintf('Skipping creation of boundary classifier for stage %d\n', stage);
    load(boundaryClassifierFilename, 'classifier');
  end

  
  fprintf('Performing merges:\n');
  for ii = 1 : consts.numImages
    if ~consts.useImages(ii)
      continue;
    end
    
    fprintf('Merging regions (Image %d/%d, stage %d)... ', ii, consts.numImages, stage);
    outFilename = sprintf(consts.boundaryInfoPostMerge, params.seg.featureSet, stage, ii);
    if exist(outFilename, 'file') && ~params.overwrite_train
      fprintf(' already exists, overwrite=false.\n');
      continue;
    end

    load(sprintf(consts.planeDataFilename, ii), 'planeData');
    load(sprintf(consts.watershedFilename, ii), 'pbAll');

    if stage == 1
      boundaryInfoFilename = sprintf(consts.watershedFilename, ii);
    else
      boundaryInfoFilename = sprintf(consts.boundaryInfoPostMerge, params.seg.featureSet, stage-1, ii);
    end
    
    load(boundaryInfoFilename, 'boundaryInfo');
    load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
    load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
    load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
    load(sprintf(consts.boundaryFeaturesFilename, params.seg.featureSet, stage, ii), 'boundaryFeatures');
    
    [~, instanceLabels] = get_labels_from_instances(boundaryInfo.imgRegions, imgObjectLabels, imgInstanceLabels);
    
    %MARKER Merging regions
    result = merge_regions(boundaryInfo, boundaryFeatures, classifier, stage, params);
    boundaryInfo = update_boundary_info(boundaryInfo, result, imgRgb, ii, stage);
    save(outFilename, 'boundaryInfo');
    fprintf(' Done.\n');
  end

  fprintf('======================================\n');
  fprintf('Finished merging regions for stage %d!\n', stage);
  fprintf('======================================\n\n');
end
