% Extracts features used for classifying image regions as one of C
% structure classes.
%
% Args:
%   ii - the image number.
%   params - the parameter struct.
%
% Returns:
%   features - RxD matrix where R is the number of regions in the image and
%              D is the dimensionality of the features.
function features = extract_region_to_structure_classes_features(ii, params)
  % fprintf('Extracting structure-class features %d\n', ii); % un-necessary output
  imgRegions = get_regions(ii, params);
  regionFeatures = extract_region_features(ii, imgRegions, params);
  siftFeatures = extract_features_sift_sc(ii, imgRegions, params);
  features = [regionFeatures, siftFeatures];
end

% Extracts a variety of region-based features.
function allRegionFeatures = extract_region_features(ii, imgRegions, params)
  Consts;
  load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
  load(sprintf(consts.imageDepthFilename, ii), 'imgDepth');
  load(sprintf(consts.planeDataFilename, ii), 'planeData');
  
  points3d = swap_YZ(planeData.points3d);
  normals = swap_YZ(planeData.normals);
 
  lowestPointInImage = min(points3d(:,3));
  
  R = max(imgRegions(:));
  
  % Create the masks only once.
  [H, W] = size(imgRegions);
  regionMasks = false(H, W, R);
  for rr = 1 : R
    regionMasks(:,:,rr) = imgRegions == rr;
  end
 
  % Finally, calculate the distance relative to the 'size' of the room. To
  % measure the size, we'll sample points in the point cloud and take as
  % the size the distance between the 2 farthest points.
  [X, Y] = meshgrid(1:50:W, 1:50:H);
  inds = sub2ind([H, W], Y(:), X(:));
  samplePoints = points3d(inds, :);
  allDists = pdist2(samplePoints, samplePoints);
  roomMaxDist = max(allDists(:));
  
  % Extract RGB features
  rgbFeatures = get_rgb_features(imgRgb, regionMasks);
  
  % Extract Surface Normal Histograms
  surfaceNormalHists = get_surf_normal_features(...
      regionMasks, points3d, normals, params);
  
  % Extract volumetric features
  volumetricFeatures = get_vol_features(...
      points3d, regionMasks);

  % Get the distance to the closest wall.
  wallPlaneIds = find(planeData.planeTypes.isVertical & ...
      planeData.planeTypes.isBoundary);
  wallPoints = cell(numel(wallPlaneIds), 1);
  
  for ii = 1 : numel(wallPlaneIds)
    wallMask = planeData.planeMap == wallPlaneIds(ii);
    
    % Take a sample of the wall points to reduce computation.
    wallPoints{ii} = points3d(wallMask, :);
    wallPoints{ii} = get_pcd_sample(wallPoints{ii}, 100, 1);
  end
  
  % Relative Depth features.
  imgMaxDepths = repmat(max(imgDepth, [], 1), [H, 1]);
  imgRelativeDepths = imgDepth ./ imgMaxDepths;
  imgDistToExtent = imgMaxDepths - imgDepth;
  
  D = 128;
  allRegionFeatures = zeros(R, D);
  
  boundingBoxes = get_bounding_boxes_2d(imgRegions);
  
  for rr = 1 : R
    sceneFeatures = get_scene_features(regionMasks(:,:,rr), points3d, ...
        imgRelativeDepths, imgDistToExtent, wallPoints, wallPlaneIds, ...
        lowestPointInImage, roomMaxDist);
    
    %%
    regionFeatures = [...
        rgbFeatures(rr,:) ...
        surfaceNormalHists(rr,:) ...
        volumetricFeatures(rr,:) ...
        boundingBoxes(rr).height ...
        boundingBoxes(rr).width ...
        sceneFeatures];

    assert(~any(isnan(regionFeatures)));
    assert(~any(isinf(regionFeatures)));
    
    allRegionFeatures(rr,:) = regionFeatures;
  end
end

% Extracts RGB features for each region. 
%
% Args:
%   imgRgb - uint8 RGB image of size HxWx3
%   imgRegions - HxW matrix, region map
%   regionIds - (optional) the subset of regions for which we want features
%               extracted.
%
% Returns:
%   features - RxD matrix where R is the number of regionIds and D is the
%              number of extracted features.
function features = get_rgb_features(imgRgb, regionMasks)
  assert(isa(imgRgb, 'uint8'));
 
  [H, W, D] = size(imgRgb);
  assert(D == 3);

  % Convert to YCbCr
  imgYCbCr = rgb2ycbcr(imgRgb);
  imgYCbCr = reshape(imgYCbCr, [H*W 3]);
  imgYCbCr = im2double(imgYCbCr);
  
  R = size(regionMasks,3);
  features = zeros(R, 36);
  
  edges = 0 : 1/10 : 1;
  edges(end) = [];
  
  for rr = 1 : R
    colorValues = imgYCbCr(regionMasks(:,:,rr), :);
    numVals = nnz(regionMasks(:,:,rr));
    features(rr, :) = [mean(colorValues) std(colorValues) ...
        histc(colorValues(:,1), edges)' / numVals ...
        histc(colorValues(:,2), edges)' / numVals...
        histc(colorValues(:,3), edges)' / numVals];
  end
end

function features = get_vol_features(points3d, regionMasks)

  R = size(regionMasks,3);
  D = 4;

  features = zeros(R, D);

  for rr = 1 : R
    regionPoints = points3d(regionMasks(:,:,rr), :);
    
    bb = get_bounding_box_3d(regionPoints);
    bbVolume = prod(2 * bb.coeffs);
    
    bbDims = sort(bb.coeffs);
    
    features(rr,1) = bbVolume;
    features(rr,2:4) = bbDims;
  end
end

function sceneFeatures = get_scene_features(regionMask, points3d, ...
    imgRelativeDepths, imgDistToExtent, wallPoints, wallPlaneIds, lowestPointInImage, roomMaxDist)

  % Take a sample of the region points to reduce computational cost.
  regionPoints = points3d(regionMask, :);
  
  if size(regionPoints,2) > 100
    regionSample = get_pcd_sample(regionPoints, 200, .9);
  else
    regionSample = regionPoints;
  end

  minDist = 10;
  for ii = 1 : numel(wallPlaneIds)
    dists = pdist2(regionSample, wallPoints{ii});
    dist = min(dists(:));
    if minDist > dist
      minDist = dist;
    end
  end
        
  meanRelativeDepth = mean(imgRelativeDepths(regionMask));
  stdDevRelativeDepth = std(imgRelativeDepths(regionMask));
  
  meanDistToExtent = mean(imgDistToExtent(regionMask));
  stdDevDistToExtent = std(imgDistToExtent(regionMask));
  
  % Max / Min Height. Use the same sampled region points to avoid using
  % outliers.
  minHeight = min(regionPoints(:,3)) - lowestPointInImage;
  maxHeight = max(regionPoints(:,3)) - lowestPointInImage;
  
  sceneFeatures = [minDist minDist / roomMaxDist ...
    meanRelativeDepth stdDevRelativeDepth ...
    meanDistToExtent stdDevDistToExtent ...
    minHeight maxHeight];
end


function siftScFeatures = extract_features_sift_sc(ii, imgRegions, params)
  Consts;
  
  % Load the Sparse Coding dictionary.
  load(sprintf(consts.siftDictionary, params.sift.patchSize, ...
      params.sift.stride, params.sift.normMethod, ...
      params.sc.K, params.sc.lambda), 'D');

  [regionIds, R] = get_region_ids(imgRegions);

  % Grab the SIFT descriptors.
  siftRgbFilename = sprintf(consts.siftRgbFilename, ...
      params.sift.patchSize, params.sift.stride, ...
      params.sift.normMethod, ii);
  
  siftDFilename = sprintf(consts.siftDepthFilename, ...
      params.sift.patchSize, params.sift.stride, ...
      params.sift.normMethod, ii);
    
  siftRgb = load(siftRgbFilename, 'features', 'coords', 'norms');
  siftD = load(siftDFilename, 'features', 'coords', 'norms');
  rgbdFeatures = [siftRgb.features siftD.features];
  coords = siftRgb.coords;

  scParamsLasso = struct();
  scParamsLasso.pos = 1;
  scParamsLasso.mode = 2;
  scParamsLasso.lambda = params.sc.lambda;
  scParamsLasso.lambda2 = 0;
  coeffs = mexLasso(rgbdFeatures', D, scParamsLasso);
  coeffs = coeffs';
    
  coeffRegions = imgRegions(sub2ind(size(imgRegions), ...
    coords(:,2), coords(:,1)));

  % Now, sum them over each region.
  siftScFeatures = zeros(R, params.sc.K);

  for rr = 1 : R
    regionCoeffs = coeffs(coeffRegions == regionIds(rr), :);
    siftScFeatures(rr,:) = sum(regionCoeffs, 1);
    siftScFeatures(rr,:) = siftScFeatures(rr,:) ./ (sum(siftScFeatures(rr,:)) + eps);
  end
end
