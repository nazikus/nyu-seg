% Performs inference by finding the optiminal binary (IP) solutions on ground truth (human
% annotated) regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

% Add Gurobi to the path.
addpath(consts.gurobiPath);
gurobi_setup;

RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));

params.regionSrc = consts.REGION_SRC_LABELS;
params.support.infMethod = consts.SUP_INF_IP;

infer_supports;
eval_supports;