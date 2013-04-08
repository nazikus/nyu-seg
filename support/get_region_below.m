% Returns the region ID of the region directly beneath the center of the
% given mask for support inference. If no region is below the given region
% but the floor IS visible in the image-plane, then the floor index is
% used.
%
% Args:
%   imgRegions - the region map, an HxW matrix whose values range from 0 to
%                R where 1..R indicate region values, 0 indicates that
%                pixel belongs to no region and R is the total number of
%                regions.
%   regionMask - the mask under which we want to search for another region.
%   canSeeFloor - whether or not the floor was found in the image plane.
%   floorRegionId - the region ID of the floor.
%
% Returns:
%   regionId - the region ID of the supporting region.
function regionId = get_region_below(imgRegions, regionMask, canSeeFloor, floorRegionId)

  DEBUG = 0;
  
  % Find the center of the object in the X axis.
  mn = find(any(regionMask,1), 1, 'first');
  mx = find(any(regionMask,1), 1, 'last');

  mid = round((mx - mn) / 2 + mn);

  % Now that we have the midpoint in X, find the first point BELOW the
  % mask.
  bottom = find(regionMask(:,mid), 1, 'last');

  % Now, find the first region below this value.
  candidates = imgRegions(bottom+1:end, mid);

  if DEBUG
    fprintf('Showing midpoint\n');
    sfigure(1);
    imagesc(regionMask);
    hold on;
    scatter(mid, bottom+1, 'g');
    hold off;
    pause;
  end
  
  if nnz(candidates) > 0
    ind = find(candidates, 1, 'first');
    regionId = candidates(ind);
    
    if DEBUG
      imgSupport = zeros(size(regionMask));
      imgSupport(regionMask) = 2;
      imgSupport(imgRegions == regionId) = 1;
      
      sfigure(1);
      imagesc(imgSupport);
      title('Showing supporting region');
      pause;
    end
  elseif canSeeFloor
    
    if DEBUG
      imgSupport = zeros(size(regionMask));
      imgSupport(regionMask) = 2;
      imgSupport(imgRegions == floorRegionId) = 1;
      
      sfigure(1);
      imagesc(imgSupport);
      title('Showing supporting region');
      pause;
    end
    
    % Supported by the observed Floor.
    regionId = floorRegionId;
  else
    % Supporting by a hidden region.
    regionId = -1;
  end
end
