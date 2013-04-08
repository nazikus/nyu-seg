% Trains the support classifier to predict the support relationship between
% two regions. Note that the input to the support classifier is NOT
% symmetric.
%
% Args:
%   params - the parameters struct. See Params.m
%   stage - the stage of segmentation (if used during segmentation) or 0 if
%           being used with the ground truth regions.
function train_support_classifier(params, stage)
  Consts;
  
  if nargin < 2
    stage = params.stage;
  end

  nn_configure_path;

  % Load the dataset.
  datasetFilename = sprintf(consts.supportDataset, params.regionSrc, stage);
  load(datasetFilename, 'trainData', 'trainLabels', 'testData', 'testLabels');

  %%
  fprintf('Normalizing data ...');
  [trainData, trainMeans] = normalize_zero_mean(trainData);
  [trainData, trainStds] = normalize_unit_var(trainData);

  testData = normalize_zero_mean(testData, trainMeans);
  testData = normalize_unit_var(testData, trainStds);
  fprintf('DONE\n');

  %%
  D = size(trainData, 2);
  C = numel(unique(trainLabels));

  nn_consts;
  nn = nn_create(D, 'support_classifier');
  nn = nn_add_layer(nn, C, SOFTMAX);

  nn.eta = 0.0001;
  nn.numUpdates = 40000;
  % nn.lambda = 0.0001;
  nn.resampleStrategy = RESAMPLE_NONE;
  % nn.resampleStrategy = RESAMPLE_RANDOM;
  nn.resampleRatio = 0;

  nn = nn_train_sgd(nn, trainData, trainLabels);

  %
  [accTrain, confMatTrain] = nn_eval(nn, trainData, trainLabels);
  [accTest, confMatTest] = nn_eval(nn, testData, testLabels);

  fprintf('AccTrain: %f\n', accTrain);
  fprintf('AccTest: %f\n', accTest);

  fprintf('\n');
  fprintf('Mean Diag train: %f\n', mean(diag(normalize_conf_mat(confMatTrain))));
  fprintf('Mean Diag test: %f\n', mean(diag(normalize_conf_mat(confMatTest))));
  fprintf('\n');

  disp(confMatTrain);
  disp(confMatTest);

  %% Save the classifier to disk.
  fprintf('Saving classifier to disk...');
  outFilename = sprintf(consts.supportClassifier, ...
      params.regionSrc, stage, nn.resampleStrategy, nn.resampleRatio);
  save(outFilename, 'nn', 'trainMeans', 'trainStds');
  fprintf('DONE.\n');

end