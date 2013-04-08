% Evaluates the segmentation pipeline using the overlap metric from Hoiem et al's Recovering Occlusion Boundaries
% from an Image.
addpath('common/');
addpath('segmentation/');
Consts;
Params;

params.seg.featureSet = consts.BFT_RGBD;

%%
load(consts.splitsPath, 'trainNdxs', 'testNdxs');

numTrain = numel(trainNdxs);
numTest = numel(testNdxs);

allWeightedScores = zeros(params.seg.numStages, 2);
allUnweightedScores = zeros(params.seg.numStages, 2);

fprintf('Evaluating Segmentations:\n');
for stage = 1 : params.seg.numStages
  
  trainInd = 0;
  testInd = 0;

  trainWeightedScores = zeros(numTrain, 1);
  trainUnweightedScores = zeros(numTrain, 1);

  testWeightedScores = zeros(numTest, 1);
  testUnweightedScores = zeros(numTest, 1);

  for ii = 1 : consts.numImages
    if ~consts.useImages(ii)
      continue;
    end
    
    fprintf('Evaluating segmentation %d/%d (stage %d)... \n', ii, consts.numImages, stage);
    % Load the segmentation regions.
    filename = sprintf(consts.boundaryInfoPostMerge, ...
        params.seg.featureSet, stage, ii);
    load(filename, 'boundaryInfo');
    imgRegionsPred = boundaryInfo.imgRegions;

    % Load the ground true regions.
    load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
    load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
    imgRegionsTrue = get_regions_from_labels(imgObjectLabels, imgInstanceLabels);
    
    [weightedScore, unweightedScore] = ...
        evaluate_segmentation(imgRegionsTrue, imgRegionsPred);

    if isin(ii, trainNdxs)
      trainInd = trainInd + 1;
      trainWeightedScores(trainInd) = weightedScore;
      trainUnweightedScores(trainInd) = unweightedScore;
    else
      testInd = testInd + 1;
      testWeightedScores(testInd) = weightedScore;
      testUnweightedScores(testInd) = unweightedScore;
    end
  end
  trainWeightedScores = trainWeightedScores(1:trainInd);
  trainUnweightedScores = trainUnweightedScores(1:trainInd);
  testWeightedScores = testWeightedScores(1:testInd);
  testUnweightedScores = testUnweightedScores(1:testInd);
  
  fprintf('\n');
  
  allWeightedScores(stage,1) = mean(trainWeightedScores);
  allWeightedScores(stage,2) = mean(testWeightedScores);
  
  allUnweightedScores(stage,1) = mean(trainUnweightedScores);
  allUnweightedScores(stage,2) = mean(testUnweightedScores);
  fprintf('done.\n');
end

%%
fprintf('======================================\n');
fprintf('              Results:\n');
fprintf('======================================\n');

for stage = 1 : params.seg.numStages
  fprintf('Stage %d:\n', stage);
  fprintf('Weighted:   Train=%2.1f  Test=%2.1f\n', ...
    100 * allWeightedScores(stage,1), ...
    100 * allWeightedScores(stage,2));
  
  fprintf('Unweighted: Train=%2.1f  Test=%2.1f\n', ...
    100 * allUnweightedScores(stage,1), ...
    100 * allUnweightedScores(stage,2));
  fprintf('\n');
end

fprintf('======================================\n');

% DEBUG
fevaluate = fopen(sprintf('%s%s_%d.txt', consts.boundaryFeaturesDir, ...
            'evaluate', consts.useNdx(end)-consts.useNdx(1)+1), 'w');

fprintf(fevaluate,['======================================\n', ...
                   '              Results:\n', ...
                   '======================================\n']);

for stage = 1 : params.seg.numStages
  fprintf(fevaluate,'Stage %d:\n', stage);
  fprintf(fevaluate,'Weighted:   Train=%2.1f  Test=%2.1f\n', ...
    100 * allWeightedScores(stage,1), ...
    100 * allWeightedScores(stage,2));
  
  fprintf(fevaluate,'Unweighted: Train=%2.1f  Test=%2.1f\n', ...
    100 * allUnweightedScores(stage,1), ...
    100 * allUnweightedScores(stage,2));
  fprintf(fevaluate,'\n');
end

fprintf(fevaluate,'======================================\n\n');
fclose(fevaluate);