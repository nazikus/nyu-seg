% Returns the ground truth support relationships for the given segmented
% regions.
%
% Args:
%   imgRegionsGt - HxW matrix, a region map with R1 regions.
%   imgRegionsSeg - HxW matrix, a region map with R2 regions.
%   supportLabelsGt - Nx3 matrix of support relationships defined on the
%                  ground truth regions.
%   imgLabels - HxW matrix, a label map for the current image.
%
% Returns:
%   supportLabelsGt - Mx3 matrix of support relationships defined on the
%                  segmented regions.
function supportLabelsSeg = get_ground_truth_support_relns(imgRegionsGt, ...
    imgRegionsSeg, supportLabelsGt, imgLabels)

  % Load the floor class ID.
  Consts;
  load(consts.datasetPath, 'namesToIds');
  floorId = namesToIds('floor');
  
  assert(size(supportLabelsGt,2) == 3);
  
  overlap = get_region_overlap(imgRegionsGt, imgRegionsSeg);  

  % Grab the ground truth support relationships for the current image. Keep
  % in mind that these are defined in terms of the ground truth regions.
  supportLabelsSeg = [];
  
  % Next, we need to remap the supports in such a way that it preserves the
  % closest supporting region.
  for ii = 1 : size(supportLabelsGt,1)
    % Find all segmented regions that map to the given ground truth
    % supportee region.
    segIdsObj = find(overlap(supportLabelsGt(ii,1),:) > .25);
    
    if isempty(segIdsObj)
      continue;
    end
    
    % First, is the region hidden or the floor?
    if supportLabelsGt(ii,2) <= 0
      for jj = 1 : numel(segIdsObj)
        segId = segIdsObj(jj);
        label = [segId supportLabelsGt(ii,2) supportLabelsGt(ii,3)];
        supportLabelsSeg = [supportLabelsSeg; label];
      end
      
      continue;
    end
    
    % Otherwise, find all segmented regions that map to the given ground
    % truth supporting region.
%     segIdsSup = find(map == supportLabelsGt(ii,2));
    segIdsSup = find(overlap(supportLabelsGt(ii,2),:) > .25);
    
    % Get an outer product of these.
    [segIdsObj, segIdsSup] = meshgrid(segIdsObj, segIdsSup);
    segIdsObj = segIdsObj(:);
    segIdsSup = segIdsSup(:);
    
    % For each segment, create a new label.
    for jj = 1 : numel(segIdsObj)
      
      segIdObj = segIdsObj(jj);
      segIdSup = segIdsSup(jj);
      
      label = [segIdObj, segIdSup, supportLabelsGt(ii,3)];
      supportLabelsSeg = [supportLabelsSeg; label];
    end
  end
  
  if isempty(supportLabelsSeg)
    return;
  end

  % Next, add any floor regions.
  regionLabels = get_labels_from_regions(imgRegionsSeg, imgLabels);
  floorRegions = find(regionLabels == floorId);
  
  for ii = 1 : numel(floorRegions)
    supportLabelsSeg = [supportLabelsSeg; [floorRegions(ii) 0 0]];
  end
  
  % Resort them in order of region ID.
  [~, inds] = sort(supportLabelsSeg(:,1), 'ascend');
  supportLabelsSeg = supportLabelsSeg(inds,:);
end
