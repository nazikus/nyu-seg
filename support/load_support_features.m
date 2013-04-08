% Loads all of the support features for a given image.
%
% Args:
%   ii - the image number.
%   params - the parameter struct.
%
% Returns:
%   features - NxD matrix of region features where n is the number of
%              region pairs in the image and D is the dimensionality of
%              each feature vector.
%   supporteeIds - Nx1 vector of region IDs.
%   supporterIds - Nx1 vector of region IDs.
function [features, supporteeIds, supporterIds] = ...
    load_support_features(ii, params, imgRegions)
  Consts;
  
  containment = load(sprintf(consts.supportFeaturesContainment, ...
      params.regionSrc, params.seg.featureSet, 0, ii), ...
        'features', 'supporteeIds', 'supporterIds');
      
  geometry = load(sprintf(consts.supportFeaturesGeometry, ...
      params.regionSrc, params.seg.featureSet, 0, ii), ...
        'features', 'supporteeIds', 'supporterIds');
  
  horz = load(sprintf(consts.supportFeaturesHorz, ...
      params.regionSrc, params.seg.featureSet, 0, ii), ...
        'features', 'supporteeIds', 'supporterIds');
      
  assert(all(containment.supporterIds == geometry.supporterIds));
  assert(all(containment.supporteeIds == geometry.supporteeIds));
  
  assert(all(containment.supporterIds == horz.supporterIds));
  assert(all(containment.supporteeIds == horz.supporteeIds));
  
  features = [containment.features geometry.features horz.features];
  supporteeIds = geometry.supporteeIds;
  supporterIds = geometry.supporterIds;

  featuresFilename = sprintf(consts.structureFeaturesFilename, ...
      params.regionSrc, params.seg.featureSet, params.stage, ii);
  load(featuresFilename, 'regionFeatures');

  D = 128 * 2 + 1;
  F = size(features,1);
  
  fullRegionFeatures = zeros(F, D);
  
  if nargin < 3
    imgRegions = get_regions(ii, params);
  end
  [regionIds, R] = get_region_ids(imgRegions);
  
  for rr = 1 : R
    N = nnz(supporteeIds == regionIds(rr));
    fullRegionFeatures(supporteeIds == regionIds(rr), 1:128) = repmat(regionFeatures(rr, 1:128), [N, 1]);
    
    N = nnz(supporterIds == regionIds(rr));
    fullRegionFeatures(supporterIds == regionIds(rr), 129:256) = repmat(regionFeatures(rr, 1:128), [N, 1]);
    
    fullRegionFeatures(supporterIds == -1, end) = 1;
  end
  
  features = [features fullRegionFeatures];
end
