% Infers the support by using some meta-class knowledge.
%
% Args:
%   imgRegions - HxW region map of the image.
%   regionIds - the region IDs for which we should infer support relations.
%   classifer - the meta-class classifier.
%   regionFeatures - features to be used by the meta-class classifier.
%
% Returns;
%   supportRelnsPred - Rx3 matrix of support relations.
%   metaClassPredictions - Rx1 vector of meta-class predictions.
function [supportRelnsPred, metaClassPredictions] = ...
  infer_supports_structure_class_rules(imgRegions, classifier, ...
    regionFeatures)

  DEBUG = 0;
  
  R = max(imgRegions(:));
  regionIds = 1 : R;
  
  regionFeatures = normalize_zero_mean(regionFeatures, classifier.trainMeans);
  regionFeatures = normalize_unit_var(regionFeatures, classifier.trainStds);
  
  metaClassProbs = nn_feed_forward(classifier.nn, regionFeatures);
  
  [~, metaClassPredictions] = max(metaClassProbs, [], 2);
  
  % Also, of the floors, figure out the most floor-y one.
  [~, floorId] = max(metaClassProbs(:,1));
  
  floorRegionId = regionIds(floorId);
  canSeeFloor = 1;
  if ~any(metaClassPredictions == 1)
    floorRegionId = -1;
    canSeeFloor = 0;
  end
  
  % Now, for each region, just choose the region directly below it as the
  % supporting region.
  supportRelnsPred = ones(R,3);
  
  for rr = 1 : R
    supportRelnsPred(rr,1) = regionIds(rr);
    
    if DEBUG
      sfigure(1);
      clf;
      imagesc(imgRegions == regionIds(rr));
      title(sprintf('Region %d', regionIds(rr)));
      pause;
    end
    
    switch metaClassPredictions(rr)
      case 1 % Floor
        supportingLabelId = 0;
        
        
        if DEBUG
          fprintf('Found floor!\n');
          sfigure(1);
          title('Its a floor!');
          pause;
        end
        
      case {2,3} % Support or Structure
        supportingLabelId = floorRegionId;
        
        if DEBUG
          imgSupport = zeros(size(imgRegions));
          imgSupport(imgRegions == regionIds(rr)) = 2;
          imgSupport(imgRegions == supportingLabelId) = 1;

          sfigure(1);
          imagesc(imgSupport);
          title('Showing supporting region');
          pause;
        end
        
      case 4 % Prop
        regionMask = imgRegions == regionIds(rr);
        supportingLabelId = get_region_below(imgRegions, regionMask, canSeeFloor, floorId);
        
        if supportingLabelId > 0 && metaClassPredictions(regionIds == supportingLabelId) == 4
          supportRelnsPred(rr,end) = 2;
        end
        
        if DEBUG
          imgSupport = zeros(size(regionMask));
          imgSupport(regionMask) = 2;
          imgSupport(imgRegions == supportingLabelId) = 1;

          sfigure(1);
          imagesc(imgSupport);
          title('Showing supporting region');
          pause;
        end
    end
    
    supportRelnsPred(rr,2) = supportingLabelId;
  end
  
  % Make sure the support type for floor/hidden is 0.
  supportRelnsPred(supportRelnsPred(:,2) <= 0,3) = 0;
end
