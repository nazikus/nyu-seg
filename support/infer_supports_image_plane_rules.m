% Infers the support relationships based on the following logic:
% 
% First, find the most likely floor-region and name it the floor.
% Next, for each other region, if it is completely surrounded on all sides
% by a single other region, then infer that its supported-from-behind by
% that region. Otherwise, label the region directly beneath it (in the
% image plane) as its support-from-below.
%
% Args:
%   imgRegions - HxW region map.
%   regionIds - the regionIds for which we should estimate support.
%   classifier - the classifier trained via run_train_floor_classifier.m
%
% Returns:
%   supportRelns - Rx3 matrix of support relations where the first column
%                  indicates supportee region IDs, the second column
%                  indicates supporter region IDs and the third column
%                  indicates support direction.
function supportRelns = infer_supports_image_plane_rules(imgRegions, ...
    classifier, regionFeatures)
  
  DEBUG = 0;
  
  R = max(imgRegions(:));
  supportRelns = zeros(R,3);
  
  canSeeFloor = 0;
  
  % First, classify each region as being floor vs no-floor.
  regionFeatures = normalize_zero_mean(regionFeatures, classifier.trainMeans);
  regionFeatures = normalize_unit_var(regionFeatures, classifier.trainStds);
  floorConfidences = nn_feed_forward(classifier.nn, regionFeatures);
  
  [floorConf, floorInd] = max(floorConfidences(:,1));
  if floorConf > .5
    canSeeFloor = 1;
  end
  
  % The area around the current region needs to be big enough in the case
  % of ground truth regions in case there is space between the selected
  % region and the neighboring regions.
  se = strel('disk', 3, 8);
  
  bb2ds = get_bounding_boxes_2d(imgRegions);
  
  for rr = 1 : R
    regionMask = imgRegions == rr;
    
    if DEBUG
      sfigure(1);
      clf;
      imagesc(regionMask);
      pause;
    end
      
    
    % First, is it the floor?
    if rr == floorInd
      
      if DEBUG
        fprintf('I think its the floor!\n');
        pause;
      end
      
      supportRelns(rr,:) = [rr 0 0];
      continue;
    end
    
    % First check if the region is completely surrounded. If so, then mark
    % it as supported horizontally.
    borderMask = imdilate(regionMask, se) & ~regionMask;
    borderRegions = unique(imgRegions(borderMask));
    borderRegions(borderRegions == 0) = [];
    
    % Make sure the bordering region is below the current region.
    
    if numel(borderRegions) == 1 && bb2ds(rr).maxY < bb2ds(borderRegions).maxY
      
      if DEBUG
        vis_selected_region(imgRegions, rr);
        fprintf('Im surrounded. Support from behind!\n');
        
        imgSupport = zeros(size(regionMask));
        imgSupport(regionMask) = 2;
        imgSupport(imgRegions == borderRegions(1)) = 1;

        sfigure(1);
        imagesc(imgSupport);
        title('Showing supporting region');
        pause;
      end
      
      supportRelns(rr,1) = rr;
      supportRelns(rr,2) = borderRegions(1);
      supportRelns(rr,3) = 2;
      continue;
    end
    
    supportingLabelId = get_region_below(imgRegions, regionMask, ...
        canSeeFloor, floorInd);

    if supportingLabelId == -1
      if canSeeFloor
        supportingLabelId = floorInd;
      else
        supportingLabelId = -1;
      end
    end
    
    supportRelns(rr,1) = rr;
    supportRelns(rr,2) = supportingLabelId;
    supportRelns(rr,3) = 1;
  end
end
