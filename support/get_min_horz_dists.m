% Returns the minimum horizontal distance between each pair of regions.
%
% Args:
%   imgRegions - the region map, an HxW matrix whose values range from 0,
%                indicating no region and R, the total number of regions in
%                the image.
%   regionIds - the regions to use.
%   points3d - Nx3 point cloud where N=H*W.
%
% Returns:
%   dists - RxR matrix where R is the number of regions.
function dists = get_min_horz_dists(imgRegions, regionIds, points3d)
  R = numel(regionIds);
  dists = zeros(R);
  
  % Sample each point cloud.
  samples = cell(R,1);
  for ii = 1 : R
    if nnz(imgRegions == regionIds(ii)) > 3
      samples{ii} = get_pcd_sample(points3d(imgRegions == regionIds(ii), :), 100, .9);
    else
      samples{ii} = points3d(imgRegions == regionIds(ii), :);
    end
  end
  
  for ii = 1 : R
    for jj = 1 : R
      dd = pdist2(samples{ii}(:,1:2), samples{jj}(:,1:2));
      dists(ii,jj) = min(dd(:));
    end
  end
end
