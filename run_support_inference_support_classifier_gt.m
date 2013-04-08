% Evaluates the accuracy of Support Baseline #3, an approach in which the
% most likely support for each region is found by using the output of a
% local classifier on ground truth regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;

params.support.infMethod = consts.SUP_INF_LCL_CLASSIFIER;

infer_supports;
eval_supports;