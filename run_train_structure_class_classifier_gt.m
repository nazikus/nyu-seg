% Trains a classifier to predict a region's Structure Class using the features extracted from the
% ground truth (manually annotated) regions.
addpath('common/');
addpath('nn/');
addpath('structure_classes/');
Consts;
Params;
nn_configure_path;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

train_structure_class_classifier(params);