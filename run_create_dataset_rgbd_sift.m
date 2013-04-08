% Creates a single (giant) dataset of all the SIFT descriptors extracted from the dataset. The
% descriptors are concatenated as per 'Indoor Scene Segmentation using a Structured Light Sensor'.
addpath('common/');
Consts;
Params;

% Setup the sample mask.
[~, sz] = get_projection_mask();

% Setup the sample mask.
sampleMask = get_sample_grid(sz(1), sz(2), ...
    params.sift.gridMargin, params.sift.stride);
  
F = nnz(sampleMask);

%% Load each descriptor.
D = 256;
N = consts.numImages;

allFeatures = zeros(N*F, D, 'single');
allNorms = zeros(N*F, 2, 'single');
allCoords = zeros(N*F, 2, 'single');
allImageNdxs = zeros(N*F, 1, 'single');

endNdx = 0;

%%
fprintf('\n');
for ii = 1 : N
  if ~consts.useImages(ii)
    continue;
  end
  fprintf('Loading sift descriptors (%d/%d)\r', ii, N);
  
  % Grab the descriptors.
  siftRgb = load(sprintf(consts.siftRgbFilename, params.sift.patchSize, ...
    params.sift.stride, params.sift.normMethod, ii), ...
      'features', 'coords', 'norms');
  siftD = load(sprintf(consts.siftDepthFilename, params.sift.patchSize, ...
    params.sift.stride, params.sift.normMethod, ii), ...
      'features', 'coords', 'norms');

  rgbdFeatures = [siftRgb.features siftD.features];
  rgbdCoords = siftRgb.coords;
  rgbdNorms = [siftRgb.norms siftD.norms];
  
  % Create the image ndxs.
  imageNdxs = ones(size(siftRgb.features, 1), 1) * ii;
  
  startNdx = endNdx + 1;
  endNdx = startNdx + size(rgbdFeatures, 1) - 1;
  
  allFeatures(startNdx:endNdx, :) = single(rgbdFeatures);
  allCoords(startNdx:endNdx, :) = single(rgbdCoords);
  allNorms(startNdx:endNdx, :) = single(rgbdNorms);
  allImageNdxs(startNdx:endNdx) = single(imageNdxs);
end
fprintf('\n');

% Truncate
allFeatures = allFeatures(1:endNdx, :);
allCoords = allCoords(1:endNdx, :);
allImageNdxs = allImageNdxs(1:endNdx);

%% Finally, save it to disk.
% outFilename = sprintf('%s\b\b\b\b_size%d.mat', ...
%     sprintf(consts.siftDataset, params.sift.patchSize, params.sift.stride, params.sift.normMethod), ...
%     length(consts.useNdx));
outFilename = sprintf(consts.siftDataset, params.sift.patchSize, params.sift.stride, params.sift.normMethod);

load(consts.splitsPath, 'trainNdxs', 'testNdxs');
allTrainIndics = isin(allImageNdxs, trainNdxs);

fprintf('Splitting into train and test...');
trainData = allFeatures(allTrainIndics, :);
testData = allFeatures(~allTrainIndics, :);

trainNorms = allNorms(allTrainIndics, :);
testNorms = allNorms(~allTrainIndics, :);
fprintf('DONE\n');
  
fprintf('Saving SIFT dataset: %s...', outFilename);
save(outFilename, 'trainData', 'testData', 'trainNorms', 'testNorms', '-v7.3');
fprintf('DONE\n');

