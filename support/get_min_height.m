% Returns the minimum height of the region after discarding potentially
% noisy depth values.
%
% Args:
%   regionMask - HxW logical matrix, the mask for the region.
%   points3d - Nx3 point cloud, where N=H*W
%
% Returns:
%   minHeight - the min height of the region (in absolute meters)
%   minInd - the index into the region.
function [minHeight, minInd] = get_min_height(regionMask, points3d)
  Z = points3d(regionMask, 3);

  % Toss out the bottom 5% of points.
  [Z, inds] = sort(Z, 'ascend');
  
  nn = round(.05 * numel(Z));
  
  minHeight = Z(nn+1);
  minInd = inds(nn+1);
end
