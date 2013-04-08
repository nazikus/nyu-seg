% Returns the features used in support classification. These include:
%
% 1. The percentage of the prop cuboid that falls over the support cuboid.
%
% Args:
%   ii - the image number.
%   imgRegions - the region map.
%   R1 - input supporteeIds.
%   R2 - input supporterIds.
%
% Returns:
%   features - NxD matrix of features where N is some number of pairs of
%              regions to consider.
function features = extract_geometry_features(ii, imgRegions, R1, R2)
  Consts;

  se = strel('disk', 10, 8);
  
  [~, R] = get_region_ids(imgRegions);
  
  load(sprintf(consts.planeDataFilename, ii), 'planeData');
  points3d = swap_YZ(planeData.points3d);
  supportSurfaces = planeData.supportSurfaces;
  
  features = [];
  supporterIds = [];
  supporteeIds = [];
  
  minHeights = zeros(R,1);
  minInds = zeros(R,1);
  
  [H, W] = size(supportSurfaces);
  masks = false(H, W, R);
  dilatedMasks = false(H, W, R);
  
  samplePoints = cell(R,1);
  numPixelsInMask = zeros(R, 1);
  
  for ii = 1 : R
    [minHeights(ii), minInds(ii)] = get_min_height(imgRegions == ii, points3d);
    masks(:,:,ii) = imgRegions == ii;
    dilatedMasks(:,:,ii) = imdilate(masks(:,:,ii), se);
    samplePoints{ii} = get_pcd_sample(points3d(masks(:,:,ii), :), 100, 1);
    numPixelsInMask(ii) = nnz(masks(:,:,ii));
  end
  
  imgHeights = reshape(points3d(:,3), size(imgRegions));
  
  minPoints3d = min(points3d(:,3));
  
  for ii = 1 : numel(R1)
    fprintf('Extracting geometry features %d/%d.\r', ii, numel(R1));
    
    srcNdx = R1(ii);
    dstNdx = R2(ii);
    
    srcMask = masks(:,:,srcNdx);
    propDilatedMask = dilatedMasks(:,:,srcNdx);
    
    propPoints = points3d(srcMask, :);
    propSamplePoints = samplePoints{srcNdx};
    
    propMinHeight = minHeights(srcNdx);
    propMinInd = minInds(srcNdx);
    
    % Be careful here because the minIntrinsicHeight may be literally EQUAL
    % to the support surface height which MAY invalidate it, so add a small
    % buffer.
    buffer = .01;
    possibleHeightMask = (imgHeights < (propMinHeight + buffer)) & ~srcMask;
    
    % By now, we've thrown away all points that are above the bottom of the
    % prop. Now throw away any point not near the bottom in the XY plane.

    % Check the vertical disparity between the object and the points directly
    % below it.
    nearbySupportMask = propDilatedMask & possibleHeightMask;
    
    numPropPixels = numPixelsInMask(srcNdx);

    supportMask = masks(:,:,dstNdx);
    supportSamplePoints = samplePoints{dstNdx};
    numSupportPixels = numPixelsInMask(dstNdx);

    supportMinHeight = minHeights(dstNdx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the two regions border one another? %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isBordering = nnz(nearbySupportMask & supportMask) > 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % How horizontal is the support region? %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numHorzPixelsInSupport = nnz(supportMask & supportSurfaces & possibleHeightMask);
    pcntHorzPixelsInSupport = numHorzPixelsInSupport / nnz(supportMask);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Whats the ratio of the size of the prop to the support? %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sizeRatio2d = numPropPixels / numSupportPixels;

    % Whats the vertical distance between the top of the supported region
    % and the bottom of the supporter region? We should NOT use the max
    % height of the support. Imagine clothing draped on a bed. If the bed
    % has bed posts, then the max(height(bed)) will be greater than any
    % clothing point. Instead, find the lowest point on the prop and find
    % the nearest support point. Measure the distance. Note that this is
    % still imperfect, in the cases when one object is supported in place
    % A but drapes over the side of the support.
    dists = pdist2(propPoints(propMinInd, :), supportSamplePoints);
    [~, closestPoint] = min(dists);

    vertDistance = propMinHeight - supportSamplePoints(closestPoint, 3);

    % Avg Height above ground. This is useful for encoding the ground.
    propHeightAboveGround = propMinHeight - minPoints3d;
    supportHeightAboveGround = supportMinHeight - minPoints3d;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % How far away are the centroids? %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    horzDistance = pdist2(mean(propSamplePoints(:,1:2)), mean(supportSamplePoints(:,1:2)));

    curFeatures = [...
      isBordering, ...
      log(numHorzPixelsInSupport+eps), ...
      pcntHorzPixelsInSupport ...
      sizeRatio2d, ...
      vertDistance, ... % The vertical distance between the min of the prop and the most likely supporting point.
      horzDistance, ...
      propHeightAboveGround, ...
      supportHeightAboveGround, ...
    ];

    features = [features; curFeatures];
  end
  
  fprintf('\n');
end
