% Evaluates the results of support inference.

% Load the train indices.
load(consts.splitsPath, 'trainNdxs');

% Load the support labels.
load(consts.supportLabels, 'supportLabels');

evalRecords = zeros(consts.numImages, 4);

confMatTrainPix = zeros(4);
confMatTestPix = zeros(4);

confMatTrainInst = zeros(4);
confMatTestInst = zeros(4);

for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  
  fprintf('Evaluating support inference %d/%d.\n', ii, consts.numImages);
  
  imgRegionsGt = get_regions(ii);
  imgRegionsSeg = get_regions(ii, params);
  
  % Load the RGB image, the structure labels and the results file for the
  % current image.
  load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
  load(sprintf(consts.structureLabelsFilename, ii), 'imgStructureLabels');
  structureLabels = get_labels_from_regions(imgRegionsSeg, imgStructureLabels);
  
  switch params.support.infMethod
    case consts.SUP_INF_IMG_PLN_RULES
      if ~exist(sprintf(consts.resultsImgFilename, params.regionSrc, ii), 'file')
        continue;
      end
      load(sprintf(consts.resultsImgFilename, params.regionSrc, ii), ...
        'supportLabelsPred');
 
    case consts.SUP_INF_STR_CLS_RULES
      if ~exist(sprintf(consts.resultsStrFilename, params.regionSrc, ii), 'file')
        continue;
      end
      load(sprintf(consts.resultsStrFilename, params.regionSrc, ii), ...
        'supportLabelsPred', 'M');
      
      imgMetaclassPred = fill_regions_with_values(imgRegionsSeg, M, 1:max(imgRegionsSeg(:)));
      validMap = imgRegionsSeg > 0 & imgStructureLabels > 0;
      [~, ~, ~, confMat] = eval_seg(imgMetaclassPred(validMap), imgStructureLabels(validMap), 4);
      
      % Accrue stats.
      hasLabel = structureLabels > 0;
      
      if isin(ii, trainNdxs)
        confMatTrainPix = confMatTrainPix + confMat;
        confMatTrainInst = confMatTrainInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      else
        confMatTestPix = confMatTestPix + confMat;
        confMatTestInst = confMatTestInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      end
      
    case consts.SUP_INF_LCL_CLASSIFIER
      if ~exist(sprintf(consts.resultsSupFilename, params.regionSrc, ii), 'file')
        continue;
      end
      load(sprintf(consts.resultsSupFilename, params.regionSrc, ii), ...
        'supportLabelsPred');
    case consts.SUP_INF_LP
      if ~exist(sprintf(consts.resultsLpFilename, params.regionSrc, ii), 'file')
        continue;
      end
      load(sprintf(consts.resultsLpFilename, params.regionSrc, ii), ...
        'supportLabelsPred', 'M');
      
      [~, M] = max(M, [], 2);
      imgMetaclassPred = fill_regions_with_values(imgRegionsSeg, M, 1:max(imgRegionsSeg(:)));
      validMap = imgRegionsSeg > 0 & imgStructureLabels > 0;
      [~, ~, ~, confMat] = eval_seg(imgMetaclassPred(validMap), imgStructureLabels(validMap), 4);
      
      % Accrue stats.
      hasLabel = structureLabels > 0;
      
      if isin(ii, trainNdxs)
        confMatTrainPix = confMatTrainPix + confMat;
        confMatTrainInst = confMatTrainInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      else
        confMatTestPix = confMatTestPix + confMat;
        confMatTestInst = confMatTestInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      end
    case consts.SUP_INF_IP
      
      if ~exist(sprintf(consts.resultsIpFilename, params.regionSrc, ii), 'file')
        continue;
      end
      load(sprintf(consts.resultsIpFilename, params.regionSrc, ii), ...
        'supportLabelsPred', 'M');
      
      [~, M] = max(M, [], 2);
      imgMetaclassPred = fill_regions_with_values(imgRegionsSeg, M, 1:max(imgRegionsSeg(:)));
      validMap = imgRegionsSeg > 0 & imgStructureLabels > 0;
      [~, ~, ~, confMat] = eval_seg(imgMetaclassPred(validMap), imgStructureLabels(validMap), 4);
      
      % Accrue stats.
      hasLabel = structureLabels > 0;
      
      if isin(ii, trainNdxs)
        confMatTrainPix = confMatTrainPix + confMat;
        confMatTrainInst = confMatTrainInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      else
        confMatTestPix = confMatTestPix + confMat;
        confMatTestInst = confMatTestInst + confusion_matrix(structureLabels(hasLabel), M(hasLabel), 4);
      end
      
  end
 

  load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
  
  supportLabelsGt = get_ground_truth_support_relns(imgRegionsGt, ...
    imgRegionsSeg, supportLabels{ii}, imgObjectLabels);
  
  % Evaluate the support predictions.
  [numMatch correctSupportTypeAgnostic, correctSupportTypeAware] = ...
    evaluate_support_predictions(supportLabelsGt, supportLabelsPred);
  
%   % Get the ground truth and inferred structure labels.
%   imgRegions = get_regions(ii, params);
%   structureLabelsGt = get_labels_from_regions(imgRegions, imgStructureLabels);
%   [~, structureLabelsLp] = max(M, [], 2);

  evalRecord = [ ...
    isin(ii, trainNdxs), ...
    numMatch, ...
    nnz(correctSupportTypeAgnostic), ...
    nnz(correctSupportTypeAware), ...
  ];
  evalRecords(ii, :) = evalRecord;
end

%%
isTrain = false(consts.numImages, 1);
isTrain(trainNdxs) = true;

% Support Type Agnostic
fprintf('\n');
fprintf('  Type       Train    Test\n');
fprintf('Agnostic:   %2.1f    %2.1f\n', ...
  100 * sum(evalRecords(isTrain, 3)) ./ sum(evalRecords(isTrain, 2)), ...
  100 * sum(evalRecords(~isTrain, 3)) ./ sum(evalRecords(~isTrain, 2)));

% Support Type Aware
fprintf('Aware:      %2.1f  %2.1f\n', ...
  100 * sum(evalRecords(isTrain, 4)) ./ sum(evalRecords(isTrain, 2)), ...
  100 * sum(evalRecords(~isTrain, 4)) ./ sum(evalRecords(~isTrain, 2)));

fprintf('\n');

%%
if params.support.infMethod == consts.SUP_INF_STR_CLS_RULES || ...
    params.support.infMethod == consts.SUP_INF_LP

  fprintf('\n');

  % Calculate accuracy.
  accTrain = sum(diag(confMatTrainPix)) / sum(confMatTrainPix(:));
  accTest = sum(diag(confMatTestPix)) / sum(confMatTestPix(:));
  fprintf('Acc Train (Pix)): %f\n', accTrain);
  fprintf('Acc Test (Pix): %f\n', accTest);

  % Calculate mean diagonal.
  fprintf('Mean diag (Train-Pix): %f\n', mean(diag(normalize_conf_mat(confMatTrainPix))));
  fprintf('Mean diag (Test-Pix): %f\n', mean(diag(normalize_conf_mat(confMatTestPix))));

  % Calculate accuracy.
  accTrain = sum(diag(confMatTrainInst)) / sum(confMatTrainInst(:));
  accTest = sum(diag(confMatTestInst)) / sum(confMatTestInst(:));
  fprintf('Acc Train (Inst)): %f\n', accTrain);
  fprintf('Acc Test (Inst): %f\n', accTest);

  % Calculate mean diagonal.
  fprintf('Mean diag (Train-Inst): %f\n', mean(diag(normalize_conf_mat(confMatTrainInst))));
  fprintf('Mean diag (Test-Inst): %f\n', mean(diag(normalize_conf_mat(confMatTestInst))));
  
  fprintf('\n');
end
