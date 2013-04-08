% Creates a dataset comprised of structure class features extracted from the segmented (bottom-up) 
% regions.
addpath('common/');
addpath('structure_classes/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;

% params.seg.featureSet = consts.BFT_RGBD_SUP_SC;
params.seg.featureSet = consts.BFT_RGBD;

create_dataset_structure_class_features(params);