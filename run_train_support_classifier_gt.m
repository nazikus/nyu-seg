addpath('common/');
addpath('nn/');
Consts;
Params;
nn_configure_path;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

% Load the dataset.
datasetFilename = sprintf(consts.supportDataset, params.regionSrc, ...
    params.seg.featureSet, params.stage);
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

nn.eta = params.support.classifier.eta;
nn.numUpdates = params.support.classifier.numUpdates;
nn.resampleStrategy = params.support.classifier.resampleStrategy;
nn.resampleRatio = params.support.classifier.resampleRatio;

nn = nn_train_sgd(nn, trainData, trainLabels);

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
outFilename = sprintf(consts.supportClassifier, params.regionSrc, ...
    params.stage, params.seg.featureSet, nn.resampleStrategy, nn.resampleRatio);
save(outFilename, 'nn', 'trainMeans', 'trainStds');
fprintf('DONE.\n');
