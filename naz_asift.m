clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 

iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
consts.matchImgId = naz_generate_match_ndxs;
consts.matchDir   = [consts.datasetDir 'asift_matched/'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'on';   
params.debug_fig = false;

% Must replace below 'consts' field when using A-SIFT

conf.imgFile = '%s/img%06d_stg%d_a';
conf.imgJuncFile      = '%s/img%06d_stg%d_b';
conf.imgSiftMatchFile = '%s/img%06d_stg%d_rt%1.1f';
conf.imgSiftHomoFile  = '%s/img%06d_stg%d_t%1.1f_rt%1.1f';

conf.sampleSize  = length(consts.useNdx); %%% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.sampleStages = [0];
conf.imgGap = 20; % size of gap between the images
conf.juncMarker = 'oy';
conf.siftMarker = 'oy';
conf.markerSize = 3;
     

% TODO to skip already computed figs

for sampleStage = conf.sampleStages

% TODO verify if classifier and corresponding hierarchical segmentation exists for given sample size
for setInd = 1:length(consts.matchImgId)

idSet = consts.matchImgId{setInd};    % set of image IDs for cuurent match group (from consts.matchImgId cell array)
imgNum = length(idSet);               % number of images in chosen set


% create subdir according to Sample size and Stage chosen
matchDir = sprintf('%ssample%d_%d/stage_%d/', consts.matchDir, floor(idSet(1)/10^5),conf.sampleSize, sampleStage);
if exist(matchDir, 'dir')~=7
     mkdir(matchDir);
end

% CONSIDERRING ONLY PAIRS OF IMAGES (BY NOW)
if imgNum < 2 || imgNum > 2; continue; end; 

% CHECK CONSISTENCY BETWEEEN SELECTED IMAGES AND USED IMAGES (const.useNdx & consts.matchImgId == naz_generate_match_ndxs)
if ~( any(consts.useNdx==idSet(1)) ||  any(consts.useNdx==idSet(2)) );  continue; end;



fprintf('\nProcessing image pairs: %d<->%d, stage %d ...', idSet(1), idSet(2), sampleStage); 

imRgb = cell(imgNum,1); % RGB images
im = cell(imgNum, 1);   % greyscale images 
x  = cell(imgNum, 1);   % junction x-coordinates
y  = cell(imgNum, 1);   % junction y-coordinates
D  = cell(imgNum, 1);   % sift descriptors

juncsIm = cell(imgNum, 1);
edgesIm = cell(imgNum, 1);

for i = 1:imgNum 
    load(sprintf(consts.imageRgbFilename,  idSet(i)), 'imgRgb');
    if sampleStage == 0
        load(sprintf(consts.watershedFilename, idSet(i)));
    else
        %load(sprintf(consts.boundaryInfoPostMerge, conf.sampleSize, params.seg.featureSet, sampleStage, idSet(1)), 'boundaryInfo');
    end
    edgesIm{i} = boundaryInfo.edges.fragments;     % edges from watershed segmentation
    juncsIm{i} = boundaryInfo.junctions.position;  % round means top-right pixel from junction point
    imgRgb = double(imgRgb)/255; % im2double
    
    %MARKER round!
    % extracting image and its junctions coordinates
    imRgb{i} = imgRgb;
    im{i} = rgb2gray(imRgb{i});
    x{i} = round(juncsIm{i}(:,1));
    y{i} = round(juncsIm{i}(:,2));
    conf.sizeX = size(imRgb{i},1);
    conf.sizeY = size(imRgb{i},2);
end
clear imgRgb boundaryInfo;
conf.imgGapStub = 0*ones(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
pairedImRgb = [imRgb{1} conf.imgGapStub imRgb{2}];

% show images
% --------------------------------------------------
filename = sprintf([conf.imgFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file')

h_img = figure('Visible', 'off'); % params.debug_visible);
imshow(pairedImRgb); axis image; axis off; title(sprintf('Images #%d #%d', idSet(1), idSet(2)));    
if params.debug; print(h_img, '-dpng', sprintf(filename, matchDir, idSet(1), sampleStage) ); end;

end

% show detected points
% --------------------------------------------------
filename = sprintf([conf.imgJuncFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file')

h_pnts = figure('Visible',params.debug_visible);
imshow(pairedImRgb); axis image; axis off; hold on; 
title(sprintf('Images #%d #%d, stage%d', idSet(1), idSet(2), sampleStage)); 
plot(juncsIm{1}(:,1),juncsIm{1}(:,2), conf.juncMarker, 'MarkerSize', conf.markerSize);
plot(juncsIm{2}(:,1)+conf.sizeY+conf.imgGap, juncsIm{2}(:,2), conf.juncMarker, 'MarkerSize', conf.markerSize);
for i = 1:imgNum
   if i==1
       shift = 0;
   else
       shift = conf.sizeY + conf.imgGap;
   end
   edges = edgesIm{i};
   for k = 1:length(edges)
        plot(edges{k}(:,1)+shift, edges{k}(:,2), 'r', 'LineWidth', 0.5);
   end
end
%if params.debug;      print(h_pnts, '-dpng', sprintf('%s/img%06d_b.png', matchDir, idSet(1)) ); end;
if params.debug;     saveas(h_pnts, filename, 'png'); end
if params.debug_fig; saveas(h_pnts, filename, 'fig'); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract descriptors (heavily blurred 21xb1 patches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf([conf.imgSiftHomoFile '.png'], matchDir, idSet(1), sampleStage, t, rt);
if ~exist(filename, 'file')

for i = 1:imgNum
   [D{i}, x{i}, y{i}] = ext_desc(im{i}, x{i}, y{i}); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute tentative matches between image 1 (a) and 2 (b) 
% by matching local features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rt        = 0.6;                % 1NN/2NN distance ratio threshold (between 0 and 1)
D2        = dist2(D{1}',D{2}'); % compute pair-wise distances between descriptors
[Y,I]     = sort(D2,2);         % sort distances
rr        = Y(:,1)./Y(:,2);     % compute D. Lowes' 1nn/2nn ratio test
inD12     = find(rr<rt);        % take only points with a 1nn/2nn ratio below 0.8
I         = I(inD12);           % select matched points
xat       = x{1}(inD12);
yat       = y{1}(inD12);
xbt       = x{2}(I);
ybt       = y{2}(I);

% show all tentative matches
h_match = figure('Visible',params.debug_visible);
imshow(pairedImRgb); axis image; axis off; hold on;
title( sprintf('Tentative correspondences: img #%d #%d, stage%d (rt=%1.1f)', idSet(1), idSet(2), sampleStage, rt) );
for i = 1:imgNum
   if i==1
       shift = 0;
   else
       shift = conf.sizeY + conf.imgGap;
   end
   edges = edgesIm{i};
   for k = 1:length(edges)
        plot(edges{k}(:,1)+shift, edges{k}(:,2), 'r', 'LineWidth', 0.5);
   end
end
plot(xat,yat,conf.siftMarker, 'MarkerSize', conf.markerSize);
plot(xbt+conf.sizeY+conf.imgGap,ybt,conf.siftMarker, 'MarkerSize', conf.markerSize);
hl = line([xat; xbt+conf.sizeY+conf.imgGap],[yat; ybt],'color','g');

if params.debug;     saveas(h_match, sprintf([conf.imgSiftMatchFile '.png'], matchDir, idSet(1), sampleStage, rt), 'png'); end
if params.debug_fig; saveas(h_match, sprintf([conf.imgSiftMatchFile '.fig'], matchDir, idSet(1), sampleStage, rt), 'fig'); end

%if length(D)<5; continue; end;


try
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robustly fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the inlier threshold (in noramlized image co-ordinates)
%t              = 0.3;
[Hab, inliers] = ransacfithomography([xat; yat], [xbt; ybt], t);

% show inliers
h_homo = figure('Visible',params.debug_visible); clf; clf;
imshow(pairedImRgb); axis image; axis off; hold on;
title(sprintf('Homography (ransac inliers): img #%d #%d, stg%d (rt=%1.1f, t=%1.1f)', idSet(1), idSet(2), sampleStage, rt, t));
for i = 1:imgNum
   if i==1
       shift = 0;
   else
       shift = conf.sizeY + conf.imgGap;
   end
   edges = edgesIm{i};
   for k = 1:length(edges)
        plot(edges{k}(:,1)+shift, edges{k}(:,2), 'r', 'LineWidth', 0.5);
   end
end
hl = line([xat(inliers); xbt(inliers)+conf.sizeY+conf.imgGap],[yat(inliers); ybt(inliers)],'color','g');
plot(xat(inliers),yat(inliers),conf.siftMarker, 'MarkerSize', conf.markerSize);
plot(xbt(inliers)+conf.sizeY+conf.imgGap,ybt(inliers), conf.siftMarker, 'MarkerSize', conf.markerSize);
if params.debug;     saveas(h_homo, filename, 'png'); end
if params.debug_fig; saveas(h_homo, filename, 'fig'); end



% ---- 
catch ME
    fprintf('Exception in homography: %s\n', ME.message);
end

else
    fprintf(' pair computed already.\n');
end

if imgNum < 3
    %pause;
    close all;
    continue
end;

% pause;
close all;
end

end

