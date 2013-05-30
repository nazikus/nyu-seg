clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 

iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
consts.matchImgId = naz_generate_match_ndxs;
consts.matchDir   = [consts.datasetDir 'asift_matched/'];

% sampleDir is used only for loading trained classifiers
consts.sampleDir   = 'sample_1449/';
consts.boundaryFeaturesDir = [consts.datasetDir  consts.boundaryDir consts.sampleDir];
consts.boundaryInfoPostMerge = [consts.boundaryFeaturesDir 's%d_info_type%d_stg%d_%06d.mat'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'on';   
params.debug_fig = false;

% Must replace below 'consts' field when using A-SIFT

conf.overwrite = true;
conf.imgFile = '%s/img%06d_stg%d_a';
conf.imgJuncFile      = '%s/img%06d_stg%d_j';
conf.imgSiftMatchFile = '%s/img%06d_stg%d_m';
conf.imgSiftMatchHori = '%s/img%06d_stg%d_mh';
conf.imgSiftMatchVert = '%s/img%06d_stg%d_mv';
conf.imgSiftHomoFile  = '%s/img%06d_stg%d_h';

conf.sampleSize  = length(consts.useNdx); %%% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.sampleStages = [5];
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



fprintf('\nProcessing image pairs: %d<->%d, stage %d ...\n', idSet(1), idSet(2), sampleStage); 

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
        load(sprintf(consts.boundaryInfoPostMerge, 1449, params.seg.featureSet, sampleStage, idSet(i)), 'boundaryInfo');
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
    
    %%% DEBUG %%%
    if i==1
        bi = boundaryInfo;
    end
    %%%%%%%%%%%%%
end
clear imgRgb boundaryInfo;
conf.imgGapStub = 0*ones(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
pairedImRgb = [imRgb{1} conf.imgGapStub imRgb{2}];

% show images
% --------------------------------------------------
filename = sprintf([conf.imgFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file') || conf.overwrite

h_img = figure('Visible', 'off'); % params.debug_visible);
imshow(pairedImRgb); axis image; axis off; title(sprintf('Images #%d #%d', idSet(1), idSet(2)));    
if params.debug; print(h_img, '-dpng', sprintf(filename, matchDir, idSet(1), sampleStage) ); end;

end

% show detected points
% --------------------------------------------------
filename = sprintf([conf.imgJuncFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file') || conf.overwrite

h_pnts = figure('Visible', params .debug_visible);
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


% applying a-sift to 2 images
% ------------- ------------------------------------
a.imIn{1} = ['imIn1.png'];
a.imIn{2} = ['imIn2.png'];
a.imgOutVert = ['imgOutVert.png'];
a.imgOutHori = ['imgOutHori.png'];
a.matchings = ['matchings.txt'];
a.keys{1} = ['keys1.txt'];
a.keys{2} = ['keys2.txt'];
a.flag_resize = 0;

tdir = [consts.datasetDir 'indoor_scene_seg_sup\ASIFT\'];
cd(tdir);
imwrite(im{1}, a.imIn{1}, 'png');
imwrite(im{2}, a.imIn{2}, 'png');
demo_ASIFT(a.imIn{1}, a.imIn{2}, a.imgOutVert, a.imgOutHori, a.matchings, a.keys{1}, a.keys{2}, a.flag_resize);
cd('..');

fid = fopen([tdir a.matchings]);
strline = regexp(fgetl(fid), ' ', 'split');
num_matches  = str2double(strline{1});
xat = zeros(num_matches, 1); yat = xat; xbt = xat; ybt = xat;
for i = 1:num_matches
    strline = regexp(fgetl(fid), ' ', 'split');
    xat(i) = str2double(strline{1});
    yat(i) = str2double(strline{3});    
    xbt(i) = str2double(strline{5});
    ybt(i) = str2double(strline{7});    
end
copyfile([tdir a.imgOutHori], ...
         sprintf([conf.imgSiftMatchHori '.png'], matchDir, idSet(1), sampleStage), 'f');
copyfile([tdir a.imgOutVert], ...
         sprintf([conf.imgSiftMatchVert '.png'], matchDir, idSet(1), sampleStage), 'f');

%TODO ~exist (its legacy from previous sift algorithm)
% % filename = sprintf([conf.imgSiftHomoFile '.png'], matchDir, idSet(1), sampleStage, t, rt);
% % if ~exist(filename, 'file')

% show all tentative matches
% ----------------------------------------------
h_match = figure('Visible',params.debug_visible);
imshow(pairedImRgb); axis image; axis off; hold on;
title( sprintf('Tentative correspondences: img #%d #%d, stage%d', idSet(1), idSet(2), sampleStage) );
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
%TODO
%hl = line([xat; xbt],[yat; ybt+conf.sizeY+conf.imgGap],'color','g');

if params.debug;     saveas(h_match, sprintf([conf.imgSiftMatchFile '.png'], matchDir, idSet(1), sampleStage), 'png'); end
if params.debug_fig; saveas(h_match, sprintf([conf.imgSiftMatchFile '.fig'], matchDir, idSet(1), sampleStage), 'fig'); end

%if length(D)<5; continue; end;
% % % 
% % % 
% % % try
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Robustly fit homography
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Specify the inlier threshold (in noramlized image co-ordinates)
% % % %t              = 0.3;
% % % [Hab, inliers] = ransacfithomography([xat; yat], [xbt; ybt], t);
% % % 
% % % % show inliers
% % % h_homo = figure('Visible',params.debug_visible); clf; clf;
% % % imshow(pairedImRgb); axis image; axis off; hold on;
% % % title(sprintf('Homography (ransac inliers): img #%d #%d, stg%d (rt=%1.1f, t=%1.1f)', idSet(1), idSet(2), sampleStage, rt, t));
% % % for i = 1:imgNum
% % %    if i==1
% % %        shift = 0;
% % %    else
% % %        shift = conf.sizeY + conf.imgGap;
% % %    end
% % %    edges = edgesIm{i};
% % %    for k = 1:length(edges)
% % %         plot(edges{k}(:,1)+shift, edges{k}(:,2), 'r', 'LineWidth', 0.5);
% % %    end
% % % end
% % % hl = line([xat(inliers); xbt(inliers)+conf.sizeY+conf.imgGap],[yat(inliers); ybt(inliers)],'color','g');
% % % plot(xat(inliers),yat(inliers),conf.siftMarker, 'MarkerSize', conf.markerSize);
% % % plot(xbt(inliers)+conf.sizeY+conf.imgGap,ybt(inliers), conf.siftMarker, 'MarkerSize', conf.markerSize);
% % % if params.debug;     saveas(h_homo, filename, 'png'); end
% % % if params.debug_fig; saveas(h_homo, filename, 'fig'); end
% % % 
% % % 
% % % 
% % % % ---- 
% % % catch ME
% % %     fprintf('Exception in homography: %s\n', ME.message);
% % % end

% % else
% %     fprintf(' pair computed already.\n');
% % end

% pause;
close all;
end

end

