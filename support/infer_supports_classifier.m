% Infers the most likely support for each region of an image by assigning
% the most likely supporting region according to a support classifier.
%
% Args:
%   imageNum - the number of the image to infer supports for.
%   params - the parameters struct. See Params.m
%
% Returns:
%   supportRelnsPred - Rx3 matrix of support relations where the first
%                      column contains supportee IDs, the second column
%                      contains supporter IDs, and the final column
%                      indicates the support direction (1 for below and 2
%                      for behind).
function supportRelnsPred = infer_supports_classifier(imageNum, params)
  Consts;

  imgRegions = get_regions(imageNum, params);
  R = max(imgRegions(:));
  
  %% Perform the support classification.
  [features, propIds, supportIds] = load_support_features(imageNum, params);
  [Dv, Dh, G] = get_local_support_probs(imageNum, features, propIds, supportIds, params);
  D = [Dv(:,1:end-1) Dh(:,1:end-1)];
  
  %%
  supportRelnsPred = ones(R,3);
  
  for rr = 1 : R
    supportRelnsPred(rr,1) = rr;
    
    if G(rr) > 0.5
      supportingLabelId = 0;
    else
      [~, supportingLabelId] = max(D(rr,:), [], 2);
      if supportingLabelId > R
        supportingLabelId = supportingLabelId - R;
        supportRelnsPred(rr,3) = 2;
      end
    end
    
    supportRelnsPred(rr,2) = supportingLabelId;
  end
end
