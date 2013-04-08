% Evaluates support prediction for a single image.
%
% Args:
%   supportRelnsTrue - the ground truth support relations, an Nx3 matrix
%                      where each row represents 
%   supportRelnsPred - the predicted labels.
%
% Returns:
%   numMatch - the number of predicted regions that are labeled.
%   correctSupportTypeAgnostic - binary vector, ignoring support type.
%   correctSupportTypeAware - binary vector, taking support type into account.
function [numMatch correctSupportTypeAgnostic, correctSupportTypeAware] = ...
      evaluate_support_predictions(supportRelnsTrue, supportRelnsPred)
  
  % Validate the input.
  assert(size(supportRelnsTrue,2) == 3);
  assert(size(supportRelnsPred,2) == 3);
  assert(min(supportRelnsTrue(:,3) >= 0));
  assert(max(supportRelnsTrue(:,3) <= 2));
  
  % First, throw out any predictions for regions for which we are missing
  % labeled support relationships.
  haveLabels = false(size(supportRelnsPred, 1), 1);
  for ii = 1 : size(supportRelnsPred)
    haveLabels(ii) = isin(supportRelnsPred(ii,1), supportRelnsTrue(:,1));
  end
  supportRelnsPred = supportRelnsPred(haveLabels, :);
  
  numMatch = size(supportRelnsPred,1);
  correctSupportTypeAgnostic = false(numMatch,1);
  correctSupportTypeAware = false(numMatch,1);
  
%   cmMetaClass = zeros(4, 4);
  
  % Next, go through each of the predictions and see if we can find them in
  % the labels (note that the true labels allow a region to be supported by
  % multiple supporting regions. Consequently, we'll count it as correct if
  % even one of these are retrieved.
  for ii = 1 : size(supportRelnsPred)
    regionId = supportRelnsPred(ii,1);
    predictedSupportId = supportRelnsPred(ii,2);
    predictedSupportType = supportRelnsPred(ii,3);
    
    % Find the true supports for the given region id.
    labelsForRegion = supportRelnsTrue(supportRelnsTrue(:,1) == regionId, :);
    
    % First, is the support region among the labels?
    labelsForRegion = labelsForRegion(labelsForRegion(:,2) == predictedSupportId,:);
    
    if ~isempty(labelsForRegion)
      correctSupportTypeAgnostic(ii) = true;
      
      % Next, is the support type correct?
      if labelsForRegion(1,3) == predictedSupportType
        correctSupportTypeAware(ii) = true;
      end
    end
  end
end
