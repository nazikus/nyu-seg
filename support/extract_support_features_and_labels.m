% Extracts support features between all pairs of regions in the image and
% saves them all to disk.
%
% Args:
%   params - the parameters struct. See Params.m
%   stage - the current stage of segmentation. This should be 0 for ground
%           truth regions and between 1-5 (inclusive) for bottom-up
%           regions. For bottom up regions, this indicates the stage from
%           which to load the regions so if you pass in stage 4 as the
%           current stage, the merged regions from stage 3 will be loaded.
%   startNdx (optional) - the image number to start with.
%   endNdx (optional) - the image number to end with.
function extract_support_features_and_labels(params, stage, ...
    startNdx, endNdx)
  error(nargchk(2,4,nargin));
  
  Consts;

  if ~exist(consts.supportFeaturesDir, 'dir')
    mkdir(consts.supportFeaturesDir);
  end
  
  if nargin < 4
    endNdx = consts.numImages;
  end
  if nargin < 3
    startNdx = 1;
  end
  
  for ii = startNdx : endNdx
    fprintf('Extracing support features %d/%d.\n', ii, consts.numImages);
    if ~consts.useImages(ii)
      continue;
    end
    
    if params.regionSrc == consts.REGION_SRC_BOTTOM_UP
      if stage == 0
        load(sprintf(consts.boundaryInfoPostMerge, ...
            params.seg.featureSet, 5, ii), 'boundaryInfo');
        imgRegions = boundaryInfo.imgRegions;
         R = max(imgRegions(:));
        [supporteeIds, supporterIds] = meshgrid(1:R, 1:R);
        supporteeIds = supporteeIds(:);
        supporterIds = supporterIds(:);
      else
        if stage <= 3
          load(sprintf(consts.boundaryInfoPostMerge, consts.BFT_RGBD, ...
            stage - 1, ii), 'boundaryInfo');
        else
          load(sprintf(consts.boundaryInfoPostMerge, params.seg.featureSet, ...
            stage - 1, ii), 'boundaryInfo');
        end
  
        % Grab the indices of the neighboring regions.
        supporteeIds = [boundaryInfo.edges.spLR(:, 1); boundaryInfo.edges.spLR(:, 2)];
        supporterIds = [boundaryInfo.edges.spLR(:, 2); boundaryInfo.edges.spLR(:, 1)];
        imgRegions = boundaryInfo.imgRegions;
      end
    else
      imgRegions = get_regions(ii);
      R = max(imgRegions(:));
      [supporteeIds, supporterIds] = meshgrid(1:R, 1:R);
      supporteeIds = supporteeIds(:);
      supporterIds = supporterIds(:);
    end

    % Extract containment features.
    features = extract_containment_features(ii, imgRegions, ...
        supporteeIds, supporterIds);
    outFilename = sprintf(consts.supportFeaturesContainment, ...
        params.regionSrc, params.seg.featureSet, stage, ii); 
    save(outFilename, 'features', 'supporteeIds', 'supporterIds');

    % Extract geometry features.
    features = extract_geometry_features(ii, imgRegions, ...
        supporteeIds, supporterIds);
    outFilename = sprintf(consts.supportFeaturesGeometry, ...
      params.regionSrc, params.seg.featureSet, stage, ii);
    save(outFilename, 'features', 'supporteeIds', 'supporterIds');

    features = extract_horz_features(ii, imgRegions, ...
        supporteeIds, supporterIds);
    outFilename = sprintf(consts.supportFeaturesHorz, ...
      params.regionSrc, params.seg.featureSet, stage, ii);
    save(outFilename, 'features', 'supporteeIds', 'supporterIds');
  end
end