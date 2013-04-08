% Infers the supporting regions, support directions and structure classes
% using the LP relaxion of the IP formulation on segmented (bottom-up) 
% regions.
addpath('common/');
addpath('nn/');
addpath('support/');
Consts;
Params;

% if consts.useGurobi
%   % Add Gurobi to the path.
%   addpath(consts.gurobiPath);
%   gurobi_setup;
% end

RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));

params.regionSrc = consts.REGION_SRC_BOTTOM_UP;
params.stage = 5;
% params.seg.featureSet = consts.BFT_RGBD_SUP_SC;
params.seg.featureSet = consts.BFT_RGBD;

params.support.infMethod = consts.SUP_INF_LP;

% infer_supports;
eval_supports;