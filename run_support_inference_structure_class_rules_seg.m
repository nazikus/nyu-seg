% Evaluates the accuracy of Support Baseline #2, a set of rules based on
% structure classes, on segmented regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;
% params.seg.featureSet = consts.BFT_RGBD_SUP_SC;
params.seg.featureSet = consts.BFT_RGBD;

params.support.infMethod = consts.SUP_INF_STR_CLS_RULES;

infer_supports;
eval_supports;