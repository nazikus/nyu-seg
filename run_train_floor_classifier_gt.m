% Trains a classifier to predict the dominant label of the classifier.
addpath('common/');
addpath('nn/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;

train_floor_classifier;