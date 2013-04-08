% Returns the local (classifier-supplied) support probabilities for a given
% image.
%
% Args:
%   imageNum - the image number.
%   features - NxD matrix of support features where N is the number of
%              pairs being considered.
%   supporteeIds - the region IDs of the regions being supported.
%   supporterIds - the region IDs of the regions doing the supporting.
%   params - the parameters struct.
%
% Returns:
%   Dv - support probabilities for support from below, an Rx(R+1) matrix
%        where row i represents the probabilitiy that region i was
%        supported from below by every other region including a hidden
%        region.
%   Dh - support probabilities for support from behind, an Rx(R+1) matrix
%        where row i represents the probability that region i was supported
%        from behind by every other region including a hidden region.
%   G - Rx1 vector for the probability that the region is the ground.
%   regionIds - the set of regionIds used.
function [Dv, Dh, G, regionIds] = get_local_support_probs(imageNum, features, ...
    supporteeIds, supporterIds, params)
  Consts;

  % Classify the support relationships
  load(sprintf(consts.supportClassifier, params.regionSrc, ...
      params.seg.featureSet, params.stage, ...
      params.support.classifier.resampleStrategy, ...
      params.support.classifier.resampleRatio), ...
    'nn', 'trainMeans', 'trainStds');
  
  features = normalize_zero_mean(features, trainMeans);
  features = normalize_unit_var(features, trainStds);
  
  preds = nn_feed_forward(nn, features);
  
  %% 
  imgRegions = get_regions(imageNum, params);
  regionIds = get_region_ids(imgRegions);
  
  regionIds = intersect([supporteeIds; supporterIds], regionIds);
  
  R = numel(regionIds);
  Rp = R + 1;

  Dv = zeros(R, Rp);
  Dh = zeros(R, Rp);
  
  % First, get all of the observed regions.
  for jj = 1 : size(preds, 1)
    if supporterIds(jj) == -1
      Dv(regionIds == supporteeIds(jj), end) = preds(jj,2);
      Dh(regionIds == supporteeIds(jj), end) = preds(jj,3);
    else
      Dv(regionIds == supporteeIds(jj), regionIds == supporterIds(jj)) = preds(jj,2);
      Dh(regionIds == supporteeIds(jj), regionIds == supporterIds(jj)) = preds(jj,3);
    end
  end
  
  %%
  G = zeros(R, 1);
  for rr = 1 : R
    G(rr) = mean(preds(regionIds(rr) == supporteeIds,4));
  end
  
  assert(~any(isnan(G)));
end
