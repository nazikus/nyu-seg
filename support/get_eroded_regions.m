% Gets the same regions eroded to avoid border issues.
%
% Args:
%   imgRegions - HxW
%   regionIds - the regions to erode.
%   sz - the number of pixels to erode.
%   minPixels - (optional) the number of minimum pixels use when eroding
%               each region. If the eroded region is left with fewer pixels
%               that this minimum, the erosion is cancelled.
function imgRegions = get_eroded_regions(imgRegions, regionIds, sz, minPixels)
  if nargin < 4
    minPixels = 10;
  end

  se = strel('disk', sz, 8);
  R = numel(regionIds);
  
  for rr = 1 : R
    regionMask = imgRegions == regionIds(rr);
    regionMask = imerode(regionMask, se);
    
    if nnz(regionMask) > minPixels
      imgRegions(imgRegions == regionIds(rr)) = 0;
      imgRegions(regionMask) = regionIds(rr);
    end
  end
end
