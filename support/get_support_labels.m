function supportRelns = get_support_labels(imgRgb, imgRegions, imgLabels, names, imageId, supportType)
  [H, W] = size(imgRegions);

  switch supportType
    case 1
      supportName = 'Vertical Support';
    case 2
      supportName = 'Horizontal Support';
    case 3
      supportName = 'Hanging Support';
    otherwise
      error('not supported');
  end
  
  supportRelns = [];
  
  candidateRegions = get_region_ids(imgRegions);
  regionIds = get_region_ids(imgRegions);
  
  while 1 == 1
    
    figure(2)
    set(gcf, 'Name', sprintf('%d: [%s]: Click prop first, then support', imageId, supportName));
    
    % Show the RGB image.
    subplot(1,2,1); imagesc(imgRgb);
    
    % Show the set of candidate support regions.
    imgCandidateRegions = get_candidate_region_img(imgRegions, candidateRegions);
    subplot(1,2,2); imagesc(imgCandidateRegions);
    
    % Grab the next candidate region.
    [x1, y1] = ginput(1);
    x1 = round(x1);
    y1 = round(y1);
    
    if x1 < 0 | y1 < 0 | x1 > W | y1 > H
      break;
%       strResponse = input('Are you done?', 's');
%       if strcmp(strResponse, 'y')
%         break;
%       end
%       
%       continue;
    end
    
    % Show the set of candidate support regions.
    imgCandidateRegions = get_candidate_region_img(imgRegions, regionIds);
    subplot(1,2,2); imagesc(imgCandidateRegions);
    
    % Grab the next supporter region.
    [x2, y2] = ginput(1);
    x2 = round(x2);
    y2 = round(y2);
    
    
    supportedRegionId = imgRegions(y1, x1);
    supporterRegionId = imgRegions(y2, x2);

    % Grab the labels. We'll need to make sure a valid region was selected.
    supportedLabel = mode(single(imgLabels(imgRegions == supportedRegionId)));
    supporterLabel = mode(single(imgLabels(imgRegions == supporterRegionId)));

    if ~isin(supporterRegionId, regionIds) 
      supportReln = [supportedRegionId -1 imageId supportType];
      supportRelns = [supportRelns; supportReln];
      fprintf('%s supported by hidden floor\n', names{supportedLabel});
    elseif supportedLabel > 0 & supporterLabel > 0
      supportReln = [supportedRegionId supporterRegionId imageId supportType];
      supportRelns = [supportRelns; supportReln];
      
      fprintf('%s supported by %s\n', names{supportedLabel}, names{supporterLabel});
    else
      fprintf('Unlabeled...\n');
    end

    candidateRegions(candidateRegions == supportedRegionId) = [];

    clear supportedLabel supporterLabel
  end
end

function imgCandidateRegions = get_candidate_region_img(imgRegions, candidateRegions)
  [H, W] = size(imgRegions);
  
  imgRegions(~isin(imgRegions, candidateRegions)) = 0;
  
  % Create the stripe mask.
  imgStripes = ones(H, W, 3);
  stride = 5;
  stripeTemplate = ones(1,W);
  stripeTemplate(1:stride:W) = 0;
  for hh = 1 : H
    shift = mod(hh-1, 5)+1;
    imgStripes(hh,:,2) = circshift(stripeTemplate', shift)';
    imgStripes(hh,:,3) = circshift(stripeTemplate', shift)';
  end
  imgStripes = reshape(imgStripes, [H*W 3]);

  se = strel('disk', 2, 8);
  mask = logical(imdilate(imgRegions, se));
  
  for ii = 1 : 3
    imgStripes(~mask,:) = 0;
  end
  
  imgCandidateRegions = imgStripes;
  
  % Now, color in each region.
  h = colormap(jet(256));
  inds = round(linspace(1, 256, numel(candidateRegions)));
  colors = h(inds, :);
  
  for ii = 1 : numel(candidateRegions)
    mask = imgRegions == candidateRegions(ii);
    imgCandidateRegions(mask, :) = repmat(colors(ii,:), [nnz(mask) 1]);
  end
  
  imgCandidateRegions = reshape(imgCandidateRegions, [H, W, 3]);
end