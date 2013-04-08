% Evaluates the accuracy of Support Baseline #1, a rule based approach
% combined with a floor classifier, on the ground truth regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

params.regionSrc = consts.REGION_SRC_LABELS;

params.support.infMethod = consts.SUP_INF_IMG_PLN_RULES;

infer_supports;
eval_supports;