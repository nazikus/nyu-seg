% Trains a classifier to predict a region's Structure Class using the features extracted from the
% segmented (bottom-up) regions.
addpath('common/');
addpath('nn/');
addpath('structure_classes/');
Consts;
Params;
% nn_configure_path;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;

% params.seg.featureSet = consts.BFT_RGBD_SUP_SC;
params.seg.featureSet = consts.BFT_RGBD;

train_structure_class_classifier(params);