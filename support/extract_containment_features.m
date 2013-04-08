% Extracts features used for classifying support relationships that relate
% to physical 'containment' in projected 2D space (onto the XY plane).
%
% Args:
%   ii - the image number.
%   imgRegions - the region map.
%   R1 - the supporteeIds
%   R2 = the supporterIds
%
% Returns:
%   features - NxD matrix of features where N is some number of pairs of
%              regions to consider.
function features = extract_containment_features(ii, imgRegions, R1, R2)
  Consts;

  load(sprintf(consts.planeDataFilename, ii), 'planeData');
  points3d = swap_YZ(planeData.points3d);
  
  supportSurfaces = planeData.supportSurfaces;
  
  [regionIds, R] = get_region_ids(imgRegions);

  features = [];
  
  [H, W] = size(supportSurfaces);
  
  minHeights = zeros(R,1);
  masks = false(H, W, R);
  samplePoints = cell(R,1);

  for ii = 1 : R
    [minHeights(ii)] = get_min_height(imgRegions == regionIds(ii), points3d);
    masks(:,:,ii) = imgRegions == ii;
    samplePoints{ii} = get_pcd_sample(points3d(masks(:,:,ii), :), 100, 1);
  end
  
  for ii = 1 : numel(R1)
    fprintf('Extracting containment features %d/%d.\r', ii, numel(R1));
    srcNdx = R1(ii);
    dstNdx = R2(ii);
    
    propMask = masks(:,:,srcNdx);

    propSamplePoints = samplePoints{srcNdx};
    
    minPropHeight = minHeights(srcNdx);
          
    supportMask = masks(:,:,dstNdx);

    % Be careful here because the minIntrinsicHeight may be literally EQUAL
    % to the support surface height which MAY invalidate it, so add a small
    % buffer.
    buffer = 0;
    possibleHeightMask = reshape(points3d(:,3) < minPropHeight + buffer, size(propMask)) & ~propMask;

    % The support points should be lower than the lowest prop point.
    supportPoints = points3d(supportMask & possibleHeightMask, :);

    % Horizontal support points must be lower than the lowest prop points.
    supportHorzPoints = points3d(supportMask & possibleHeightMask & supportSurfaces, :);

    if size(supportPoints, 1) <= 3
      percOfPropAnySupported = 0;
    else
      percOfPropAnySupported = get_pcnt_contained_2d(propSamplePoints(:,1:2), supportPoints(:,1:2));
    end

    if percOfPropAnySupported > 0 && size(supportHorzPoints,1) > 10
      percOfPropHorzSupported = get_pcnt_contained_2d(propSamplePoints(:,1:2), supportHorzPoints(:,1:2));
    else
      percOfPropHorzSupported = 0;
    end

    curFeatures = [percOfPropAnySupported percOfPropHorzSupported];

    features = [features; curFeatures];
  end
end
