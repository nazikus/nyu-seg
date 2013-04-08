% Extract SIFT features from each frame's RGB and Depth image.
addpath('common/');
Consts;
Params;

OVERWRITE = params.overwrite_sift;

[~, sz] = get_projection_mask();

% Setup the sample mask.
sampleMask = get_sample_grid(sz(1), sz(2), ...
    params.sift.gridMargin, params.sift.stride);
[Y, X] = ind2sub(size(sampleMask), find(sampleMask));
coords = [X(:) Y(:)];

if ~exist(consts.siftDir, 'dir')
  mkdir(consts.siftDir);
end

%%
for ii = 1 : consts.numImages
  if ~consts.useImages(ii)
    continue;
  end
  fprintf('Extracting SIFT descriptors %d/%d... ', ii, consts.numImages);
  
  % Extract from RGB.
  rgbFilename = sprintf(consts.siftRgbFilename, params.sift.patchSize, ...
      params.sift.stride, params.sift.normMethod, ii);
  if ~exist(rgbFilename, 'file') || OVERWRITE
    load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
    imgGray = rgb2gray(im2double(imgRgb));
    [features, norms] = extract_sift(imgGray, coords, params.sift);
    save(rgbFilename, 'features', 'coords', 'norms');
  else
      fprintf('rgb sifts already exists, overwrite=false');
  end
  
  depthFilename = sprintf(consts.siftDepthFilename, params.sift.patchSize, ...
      params.sift.stride, params.sift.normMethod, ii);
  if ~exist(depthFilename, 'file') || OVERWRITE
    load(sprintf(consts.imageDepthFilename, ii), 'imgDepth');

    % Make the depth relative.
    imgDepth = imgDepth - min(imgDepth(:));
    imgDepth = imgDepth ./ max(imgDepth(:));

    [features, norms] = extract_sift(imgDepth, coords, params.sift);
    save(depthFilename, 'features', 'coords', 'norms');
  else
     fprintf(', depth sifts already exists, overwrite=false');
  end
  fprintf('.\n');
end

