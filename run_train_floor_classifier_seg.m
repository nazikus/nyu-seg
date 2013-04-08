% Trains a classifier to predict the dominant label of the classifier.
addpath('common/');
addpath('nn/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;

% params.seg.featureSet = consts.BFT_RGBD_SUP_SC;
params.seg.featureSet = consts.BFT_RGBD;

train_floor_classifier;