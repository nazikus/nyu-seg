% Creates a dataset comprised of structure class features extracted from the ground truth (manually
% annotated) regions.
addpath('common/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

create_dataset_structure_class_features(params);