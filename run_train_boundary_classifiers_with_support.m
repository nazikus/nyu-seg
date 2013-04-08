% Trains several stages of boundary classifiers.
addpath('common/');
addpath('nn/');
addpath('segmentation/');
addpath('structure_classes/');
addpath('support/');
addpath(genpath('iccv07Final'));

Consts;
addpath(consts.spamsPath);

Params;
params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.seg.featureSet = consts.BFT_RGBD_SUP;

% Load the train/test split.
load(consts.splitsPath, 'trainNdxs');

if ~exist(consts.boundaryFeaturesDir, 'dir')
  mkdir(consts.boundaryFeaturesDir);
end

OVERWRITE = false;

%%
for stage = 3 : params.seg.numStages
  params.stage = stage;
  
  if stage >= 3 && ...
      (params.seg.featureSet == consts.BFT_RGBD_SUP || ...
       params.seg.featureSet == consts.BFT_RGBD_SUP_SC)
    
    dummyParams = params;
    dummyParams.stage = stage - 1;           
     
    for ii = 1 : consts.numImages
      % Extract structure class features from the regions at the previous
      % stage.
      
      outFilename = sprintf(consts.structureFeaturesFilename, ...
        dummyParams.regionSrc, dummyParams.seg.featureSet, dummyParams.stage, ii);
      if exist(outFilename, 'file') && ~OVERWRITE
        continue;
      end

      regionFeatures = extract_region_to_structure_classes_features(ii, dummyParams);
      save(outFilename, 'regionFeatures');
    end
    
    %%
    create_dataset_structure_class_features(dummyParams);
    train_structure_class_classifier(dummyParams);
   
    %%
    extract_support_features_and_labels(dummyParams, stage);
    
    %%
    create_dataset_support_features_for_seg(dummyParams, stage);
    train_support_classifier(dummyParams, stage);
  end
  
  %%
  extract_boundary_classifier_features_and_labels(stage, params);

  %% Create the boundary-classification dataset.
  datasetFilename = sprintf(consts.boundaryFeaturesDataset, ...
      params.seg.featureSet, stage);

  if ~exist(datasetFilename, 'file') || params.overwrite
    [trainData, testData, trainLabels, testLabels] = ...
        create_boundary_classifier_dataset(stage, trainNdxs, params);
    fprintf('Saving dataset...');
    save(datasetFilename, 'trainData', 'trainLabels', ...
        'testData', 'testLabels', '-v7.3');
    fprintf('DONE\n');
  else
    fprintf('Loading the boundary-classification dataset.\n');
    load(datasetFilename, 'trainData', 'trainLabels', ...
      'testData', 'testLabels');
  end

  %% Train the boundary classifier.
  boundaryClassifierFilename = ...
      sprintf(consts.boundaryClassifierFilename, params.seg.featureSet, stage);

  if ~exist(boundaryClassifierFilename, 'file') || params.overwrite
    classifier = train_boundary_classifier_dt(stage, trainData, trainLabels, ...
        testData, testLabels, params);
    save(boundaryClassifierFilename, 'classifier');
  else
    fprintf('Skipping creation of boundary classifier for stage %d\n', stage);
    load(boundaryClassifierFilename, 'classifier');
  end

  %%
  fprintf('\n');
  for ii = 1 : consts.numImages
    fprintf('Performing merge for image %d/%d (stage %d).\n', ...
        ii, consts.numImages, stage);

    if ~consts.useImages(ii)
      continue;
    end
    
    outFilename = sprintf(consts.boundaryInfoPostMerge, ...
          params.seg.featureSet, stage, ii);
    if exist(outFilename, 'file') && ~params.overwrite
      continue;
    end

    load(sprintf(consts.planeDataFilename, ii), 'planeData');
    load(sprintf(consts.watershedFilename, ii), 'pbAll');

    if stage == 1
      boundaryInfoFilename = sprintf(consts.watershedFilename, ii);
    elseif stage <= 3 && ...
        (params.seg.featureSet == consts.BFT_RGBD_SUP || ...
         params.seg.featureSet == consts.BFT_RGBD_SUP_SC)
      boundaryInfoFilename = sprintf(consts.boundaryInfoPostMerge, ...
          consts.BFT_RGBD, stage-1, ii);
    else
      boundaryInfoFilename = sprintf(consts.boundaryInfoPostMerge, ...
          params.seg.featureSet, stage-1, ii);
    end
    
    load(boundaryInfoFilename, 'boundaryInfo');
    load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
    load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
    load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
    load(sprintf(consts.boundaryFeaturesFilename, ...
        params.seg.featureSet, stage, ii), 'boundaryFeatures');
    
    [~, instanceLabels] = get_labels_from_instances(boundaryInfo.imgRegions, ...
        imgObjectLabels, imgInstanceLabels);
    
    result = merge_regions(boundaryInfo, boundaryFeatures, ...
        classifier, stage, params);
    boundaryInfo = update_boundary_info(boundaryInfo, result, imgRgb);
    save(outFilename, 'boundaryInfo');
  end
  fprintf('\n');
end
