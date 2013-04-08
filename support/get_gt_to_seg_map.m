% Gets the mapping from the ground truth regions to the segmented (bottom
% up) regions.
%
% Args:
%   imgRegionsGt - ground truth region map, a HxW matrix where H and W are
%                  the width and height of the region map, respectively.
%   imgRegionsSeg - segmentation region map, a HxW matrix where H and W are
%                  the width and height of the region map, respectively.
%
% Returns:
function map = get_gt_to_seg_map(imgRegionsGt, imgRegionsSeg)
  
  R = max(imgRegionsSeg(:));
  
  overlap = get_region_overlap(imgRegionsGt, imgRegionsSeg);
  
  map = zeros(R, 1);
  
  % For each Ground Truth region, find the closest matching Segmented
  % region.
  for rr = 1 : R
    [~, bestRegionIdPred] = max(overlap(:,rr));
    
    if overlap(bestRegionIdPred, rr) >= .25
      map(rr) = bestRegionIdPred;
    end    
  end
  assert(max(map) <= max(imgRegionsGt(:)));
end