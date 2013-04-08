% Extracts support features for every image in the dataset using ground
% truth regions.
addpath('common/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;
params.stage = 0;
params.seg.featureSet = 0;

% Whether or not to overwrite the support file if it already exists.
OVERWRITE = false;

extract_support_features_and_labels(params, 0);
