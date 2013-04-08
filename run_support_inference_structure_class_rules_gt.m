% Evaluates the accuracy of Support Baseline #2, a set of rules based on
% structure classes, on ground truth regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;

params.support.infMethod = consts.SUP_INF_STR_CLS_RULES;

infer_supports;
eval_supports;