% Gets the initial minimum heights of each region.
%
% Args:
%   imgRegions - HxW image.
%   regionIds - Rx1 vector
%   points3d - Nx3 matrix where N=H*W
%
% Returns:
%   minHeights - Rx1 vector.
%   maxHeights - Rx1 vector.
function [minHeights, maxHeights] = get_initial_heights(imgRegions, regionIds, points3d)
  imgRegions = get_eroded_regions(imgRegions, regionIds, 6, 10);
  
  % Toss out top/bottom 5%?
  R = numel(regionIds);
  minHeights = zeros(R, 1);
  maxHeights = zeros(R, 1);
  
  for rr = 1 : R
    regionMask = imgRegions == regionIds(rr);
    hh = points3d(regionMask,3);
    minHeights(rr) = min(hh);
    maxHeights(rr) = max(hh);
  end
end
