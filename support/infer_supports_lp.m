% Infers the support relations and structure-classes using an LP relaxion
% of the IP formuation for support inference.
%
% Args:
%   imageNum - the image number for which we are inferring support.
%   params - the parameters struct. See Params.m
%
% Returns:
%   supportRelnsPred - the predicted support relations.
%   S - the S-matrix, a 0-1 representation of support, an Rx(2(R+1)+1)
%       matrix.
%   M - the support matrix, an Rx4 matrix whose rows sum to 1 and whose
%       columns indicate the assignment to each of the structure classes.
%   E - the energy of the solution.
function [supportRelnsPred, S, M, E] = infer_supports_lp(imageNum, params)
  Consts;

  imgRegions = get_regions(imageNum, params);
  
  % Load the containment features.
  [features, propIds, supportIds] = load_support_features(imageNum, params, imgRegions);
  [Dv, Dh, G, regionIds] = get_local_support_probs(imageNum, features, propIds, supportIds, params);
  D  = [Dv Dh G];
  
  % Add cost for 'Hidden Floor'.
  R = numel(regionIds);
  Rp = R + 1;
  
  T = get_structure_class_distributions(imageNum, params);
  
  % Get the min heights and add one for the missing floor.
  load(sprintf(consts.planeDataFilename, imageNum), 'planeData');
  points3d = swap_YZ(planeData.points3d);
  [Hb, Ht] = get_initial_heights(imgRegions, regionIds, points3d);
  Hb = [Hb; min(Hb)];
  Ht = [Ht; min(Hb)];
  
  % Get the min horz distances and add one for the missing floor.
  V = get_min_horz_dists(imgRegions, regionIds, points3d);
  V = [V zeros(R,1)];
  
  [TC_V, TC_H] = get_support_type_trans_costs(params);

  unaryCostSupport = double(-log(D + eps));
  unaryCostMetaclass = double(-log(T + eps));
  
  switch params.regionSrc
    case consts.REGION_SRC_LABELS
      coeffs = params.support.coeffsGt;
    case consts.REGION_SRC_BOTTOM_UP
      coeffs = params.support.coeffsSeg;
  end
  
  unaryCostSupport(:,Rp) = coeffs.rho2 / R;
  
  [S, M, E] = minimize_support_energy(unaryCostSupport, unaryCostMetaclass, ...
      TC_V, TC_H, Hb, Ht, V, params);
  
  supportRelnsPred = get_support_relns_from_S(S, regionIds);
end

function supportRelnsPred = get_support_relns_from_S(S, regionIds)
  R = numel(regionIds);

  [~, ss] = max(S, [], 2);
  
  supportRelnsPred = zeros(R,3);
  supportRelnsPred(:,1) = regionIds;
  
  for ii = 1 : R
    if ss(ii) == 2*(R+1)+1
      supportRelnsPred(ii,2) = 0;
      continue;
    elseif ss(ii) == 2*(R+1) || ...
        ss(ii) == R+1
      supportRelnsPred(ii,2) = -1;
      continue;
    end
      
    if ss(ii) > R + 1
      supportRelnsPred(ii,3) = 2;
      ss(ii) = ss(ii) - (R+1);
    else
      supportRelnsPred(ii,3) = 1;
    end
    
    supportRelnsPred(ii,2) = regionIds(ss(ii));
  end
end

function T = get_structure_class_distributions(ii, params)
  Consts;

  % Produce the distribution of Region->SupportType
  load(sprintf(consts.structureClassifier, params.regionSrc, ...
      params.seg.featureSet, params.stage), ...
    'nn', 'trainMeans', 'trainStds');
  
  load(sprintf(consts.structureFeaturesFilename, params.regionSrc, ...
      params.seg.featureSet, params.stage, ii), 'regionFeatures');
  
  regionFeatures = normalize_zero_mean(regionFeatures, trainMeans);
  regionFeatures = normalize_unit_var(regionFeatures, trainStds);
  T = nn_feed_forward(nn, regionFeatures);
end
