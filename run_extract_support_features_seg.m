% Extracts support features for every image in the dataset using ground
% truth regions.
addpath('common/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;
params.seg.featureSet = consts.BFT_RGBD_SUP_SC;

% Whether or not to overwrite the support file if it already exists.
OVERWRITE = false;

extract_support_features_and_labels(params, 5);
