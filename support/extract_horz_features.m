% Extracts features useful for determining support-from-behind.
%
% Args:
%   ii - the image number.
%   imgRegions - the region map.
%   R1 - the source regions.
%   R2 - the destination regions.
%
% Returns:
%   features - NxD matrix of features where N is some number of pairs of
%              regions to consider.
function features = extract_horz_features(ii, imgRegions, R1, R2)
  Consts;
  
  load(sprintf(consts.planeDataFilename, ii), 'planeData');
  points3d = swap_YZ(planeData.points3d);
  normals = swap_YZ(planeData.normals);
  
  [H, W] = size(planeData.supportSurfaces);
  
  normals = flip_normals_towards_interior(points3d, normals);
  
  [~, inclinations, azimuths] = cart2sphere(normals(:,1), normals(:,2), normals(:,3));
  inclinations = reshape(inclinations, size(imgRegions));
  azimuths = reshape(azimuths, size(imgRegions));
  azimuths(azimuths > pi) = azimuths(azimuths > pi) - pi;
  
  distToHorz = abs(pi/2 - inclinations);

  R = max(imgRegions(:));

  masks = false(H,W,R);
  dilatedMasks = false(H,W,R);
  
  numPixelsInRegion = zeros(R, 1);
  
  pcntBehindCentroid = zeros(R,1);
  candidateRegionIds = cell(R, 1);
  numCandidateRegionIds = zeros(R, 1);
  numVertPixelsPerRegion = zeros(R, 1);
  pcntVertPixelsPerRegion = zeros(R, 1);
  
  se = strel('disk', 9, 8);
  
  numAzBins = 8;
  azimuthBins = linspace(0, pi, numAzBins);
  azimuthHists = zeros(R, numAzBins);
  
  for ii = 1 : R
    mask = imgRegions == ii;
    masks(:,:,ii) = mask;
    dilatedMasks(:,:,ii) = imdilate(mask, se);
    
    numPixelsInRegion(ii) = nnz(mask);
    
    % Calculate the 'hanging' feature.
    regionPoints = points3d(mask,:);
    bb2d = get_bounding_box_2d(mask);
    
    regionMask2 = dilatedMasks(:,:,ii) & ~mask;
    regionMask2(1:bb2d.maxY,:) = false;
    
    belowPoints = points3d(regionMask2,:);
    
    % Now, check to see if these points are farther than the centroid of
    % the bounding 3d volume.
    bb3d = get_bounding_box_3d(regionPoints);
    
    dists = sqrt(sum(belowPoints(:,1:2).^2,2));
    distCentroid = sqrt(sum(bb3d.centroid(1:2).^2));
    numFarther = nnz(dists > distCentroid);
    
    if isempty(dists)
      pcntBehindCentroid(ii) = 0;
    else
      pcntBehindCentroid(ii) = numFarther / numel(dists);
    end
    
    % For horizontal supports, look around it.
    nearbySupportMask = dilatedMasks(:,:,ii) & ~mask;
    
    % Now, for each of the regions in the possible height mask, collect
    % features that may indicate vertical support.
    candidateRegionIds{ii} = unique(imgRegions(nearbySupportMask));
    candidateRegionIds{ii}(candidateRegionIds{ii} == 0) = [];
    
    numCandidateRegionIds(ii) = numel(candidateRegionIds{ii});
    
    numVertPixelsPerRegion(ii) = nnz(distToHorz < .4 & mask);
    pcntVertPixelsPerRegion(ii) = numVertPixelsPerRegion(ii) ./ numPixelsInRegion(ii);
    
    azimuthHists(ii,:) = histc(azimuths(mask), azimuthBins);
    azimuthHists(ii,:) = azimuthHists(ii,:) ./ sum(azimuthHists(ii,:));
  end

  features = [];
  
  for ii = 1 : numel(R1)
    fprintf('Extracting horz features %d/%d.\r', ii, numel(R1));
    
    srcNdx = R1(ii);
    dstNdx = R2(ii);
    
    srcNumVertPixels = numVertPixelsPerRegion(srcNdx);
    srcPcntVertPixels = pcntVertPixelsPerRegion(srcNdx);

    dstNumVertPixels = numVertPixelsPerRegion(dstNdx);
    dstPcntVertPixels = pcntVertPixelsPerRegion(dstNdx);

    isCand = isin(dstNdx, candidateRegionIds{srcNdx});
    histDiff = chi_squ_dist(azimuthHists(srcNdx,:), azimuthHists(dstNdx,:));

    curFeatures = [...
      log(srcNumVertPixels+eps), ... % Verticallity of source region.
      srcPcntVertPixels, ...
      log(dstNumVertPixels+eps), ... % Verticallity of target region.
      dstPcntVertPixels, ...
      numCandidateRegionIds(srcNdx) ...
      isCand ...
      pcntBehindCentroid(srcNdx) ...
      histDiff ...
    ];

    features = [features; curFeatures];
  end
  
  fprintf('\n');
end

function dist = chi_squ_dist(x, y)
  aa = (x - y).^2;
  bb = x+y;
  aa(bb > 0) = aa(bb > 0) ./ bb(bb > 0);
  dist = sum(aa) / 2;
end
