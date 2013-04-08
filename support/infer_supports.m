% Script that runs the support inference.
nn_configure_path;

% Load the support labels.
load(consts.supportLabels, 'supportLabels');

% Load the train/test split.
load(consts.splitsPath, 'trainNdxs', 'testNdxs');

%%
for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  
  fprintf('Running inferring support %d/%d.\n', ii, consts.numImages);
  
  % Load the ground truth regions.
  imgRegionsGt = get_regions(ii);
  
  % Load the predicted regions.
  imgRegionsSeg = get_regions(ii, params);

  % Load the ground truth support relationships.
  load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
  
  supportRelnsTrue = get_ground_truth_support_relns(imgRegionsGt, ...
    imgRegionsSeg, supportLabels{ii}, imgObjectLabels);
  
  if isempty(supportRelnsTrue)
    continue;
  end
  
  %% Next, predict the support relationships.
  switch params.support.infMethod
    case consts.SUP_INF_IMG_PLN_RULES
      % Load the Floor classifier.
      classifier = load(sprintf(consts.floorClassifier, ...
        params.regionSrc, params.stage), 'nn', 'trainMeans', 'trainStds');
      
      % Load the Structure-class Features.
      outFilename = sprintf(consts.structureFeaturesFilename, ...
          params.regionSrc, params.seg.featureSet, params.stage, ii);
      load(outFilename, 'regionFeatures');

      % Predict the supports.
      supportLabelsPred = infer_supports_image_plane_rules(imgRegionsSeg, ...
          classifier, regionFeatures);
        
      save(sprintf(consts.resultsImgFilename, params.regionSrc, ii), ...
          'supportLabelsPred');
    case consts.SUP_INF_STR_CLS_RULES
      % Load the Structure-Class classifier.
      classifier = load(sprintf(consts.structureClassifier, ...
        params.regionSrc, params.seg.featureSet, params.stage), ...
          'nn', 'trainMeans', 'trainStds');
        
      % Load the Structure-class Features.
      outFilename = sprintf(consts.structureFeaturesFilename, ...
          params.regionSrc, params.seg.featureSet, params.stage, ii);
      load(outFilename, 'regionFeatures');
      
      % Predict the supports.
      [supportLabelsPred, M] = infer_supports_structure_class_rules(...
          imgRegionsSeg, classifier, regionFeatures);
      save(sprintf(consts.resultsStrFilename, params.regionSrc, ii), ...
          'supportLabelsPred', 'M');
    case consts.SUP_INF_LCL_CLASSIFIER
      supportLabelsPred = infer_supports_classifier(ii, params);
      save(sprintf(consts.resultsSupFilename, params.regionSrc, ii), ...
          'supportLabelsPred');
    case consts.SUP_INF_LP
      % Next, predict the support relationships.
      [supportLabelsPred, S, M, E_LP] = infer_supports_lp(ii, params);
      save(sprintf(consts.resultsLpFilename, params.regionSrc, ii), ...
        'supportLabelsPred', 'S', 'M', 'E_LP');
    case consts.SUP_INF_IP
      % Next, predict the support relationships.
      [supportLabelsPred, S, M, E_IP] = infer_supports_lp(ii, params);
      save(sprintf(consts.resultsIpFilename, params.regionSrc, ii), ...
        'supportLabelsPred', 'S', 'M', 'E_IP');
  end
end
