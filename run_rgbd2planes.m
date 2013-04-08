% Finds major scene surfaces and rotates the entire scene such that the
% floor plane's surface normal points directly up.
addpath('common');
addpath(genpath('iccv07Final'));
addpath(genpath('graph_cuts'));
addpath('surfaces');

Consts; Params;
% Whether or not to overwrite the planeData files if they're already found
% on disk.
OVERWRITE = params.overwrite_planes;

if ~exist(consts.planeDataDir, 'dir')
  mkdir(consts.planeDataDir);
end
global ii_; 
%
for ii_ = 1 : consts.numImages
%  fprintf('Extracting plane data (%d/%d)... ', ii_, consts.numImages);

  if ~consts.useImages(ii_)
%    fprintf('skipping.\n');
    continue;
  end
  
  fprintf('Extracting plane data (%d/%d)... ', ii_, consts.numImages);
  outFilename = sprintf(consts.planeDataFilename, ii_);
  if exist(outFilename, 'file') && ~OVERWRITE
    fprintf('skipping (exists), OVERWRITE=false.\n');
    continue
  end
  
  load(sprintf(consts.imageRgbFilename, ii_), 'imgRgb');
  load(sprintf(consts.imageDepthFilename, ii_), 'imgDepthOrig');
  load(sprintf(consts.imageDepthRawFilename, ii_), 'imgDepthRawOrig');
  load(sprintf(consts.surfaceNormalData, ii_), 'imgNormals', 'normalConf');
  
  planeData = rgbd2planes(imgRgb, imgDepthOrig, imgDepthRawOrig, ...
      imgNormals, normalConf);
    
  save(outFilename, 'planeData');
  fprintf('done!\n');
end
clearvars -global ii_;

fprintf('\n\n ================================================= \n\n');
