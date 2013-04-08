% Trains the structure class classifier and saves it to disk.
%
% Args:
%   params - the parameter struct to guide training.
function train_structure_class_classifier(params)
  Consts;
  
  OVERWRITE = params.overwrite_trainstruct;

  % Load the train test split.
  load(consts.splitsPath, 'trainNdxs');

  %% Load the train/test sets.
  fprintf('Loading dataset...');
  datasetFilename = sprintf(consts.structureFeaturesDataset, params.regionSrc, params.seg.featureSet, params.stage);
  load(datasetFilename, 'trainData', 'trainLabels', 'testData', 'testLabels');
  fprintf('DONE.\n');

  C = numel(unique(trainLabels));

  %%
  fprintf('Normalizing data ...');
  [trainData, trainMeans] = normalize_zero_mean(trainData);
  [trainData, trainStds] = normalize_unit_var(trainData);

  testData = normalize_zero_mean(testData, trainMeans);
  testData = normalize_unit_var(testData, trainStds);
  fprintf('DONE\n');

  D = size(trainData, 2);

  %%
  RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));
  lambdas = .0001;

  nn_consts;
  nn = nn_create(D, 'region_classifier');
  nn = nn_add_layer(nn, C, SOFTMAX);

  nn.eta = 0.0005;
  nn.numUpdates = 5000;
  nn.lambda = lambdas;


  outFilename = sprintf(consts.structureClassifier, params.regionSrc, params.seg.featureSet, params.stage);
  % outFilename = sprintf('%s\b\b\b\b_size%d.mat', sprintf(consts.structureClassifier, params.regionSrc, params.seg.featureSet, params.stage), length(consts.useNdx));
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('\n --- Skipping training --- already trained (%s).\n', outFilename);    
  else
      %DEBUG
      global debug_stream;
      if params.debug
        debug_stream = fopen([outFilename '.txt'], 'w'); 
      else
        debug_stream = 1; % 1 == stdout
      end
      %=====

      nn = nn_train_sgd(nn, trainData, trainLabels);

      % Evaluate
      [accTrain, cmTrain, ranksTrain] = nn_eval(nn, trainData, trainLabels);
      [accTest, cmTest, ranksTest] = nn_eval(nn, testData, testLabels);
      fprintf('Acc Train: %f\n', accTrain);
      fprintf('Acc Test: %f\n', accTest);

      fprintf('Mean diag (Train): %f\n', mean(diag(normalize_conf_mat(cmTrain))));
      fprintf('Mean diag (Test): %f\n', mean(diag(normalize_conf_mat(cmTest))));

      %DEBUG
      fprintf(debug_stream, '\n\n============================\n\n');
      fprintf(debug_stream, 'Acc Train: %f\n', accTrain);
      fprintf(debug_stream, 'Acc Test: %f\n', accTest);

      fprintf(debug_stream, 'Mean diag (Train): %f\n', mean(diag(normalize_conf_mat(cmTrain))));
      fprintf(debug_stream, 'Mean diag (Test): %f\n', mean(diag(normalize_conf_mat(cmTest))));
      
      if params.debug
        fclose(debug_stream);
      end
      clearvars -global debug_stream;
      %=====

      % Save the results and the normalization variables to disk.
      fprintf('Saving classifier to file %s...', outFilename);
      save(outFilename, 'nn', 'cmTrain', 'cmTest', 'trainMeans', 'trainStds', '-v7.3');
      fprintf('DONE.\n');
  end

end