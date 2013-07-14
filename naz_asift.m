clear; sclose all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 
warning off all;

iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
%consts.matchImgId = naz_generate_match_ndxs;
consts.matchDir   = [consts.datasetDir 'asift_matched/'];

% sampleDir is used only for loading trained classifiers
consts.sampleDir   = 'sample_1449/'; % meaning that 1449 classifier will be always used
consts.boundaryFeaturesDir = [consts.datasetDir  consts.boundaryDir consts.sampleDir];
consts.boundaryInfoPostMerge = [consts.boundaryFeaturesDir 's%d_info_type%d_stg%d_%06d.mat'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'off';   
params.debug_fig = false;

% Must replace below 'consts' field when using A-SIFT
conf.overwrite_image = true;
conf.overwrite_asift = true;
conf.imgFile = '%s/img%06d_stg%d_a';
conf.imgJuncFile       = '%s/img%06d_stg%d_j';         % edges and junctions
conf.imgSiftMatchPoint = '%s/img%06d_stg%d_l';         % asift match points
conf.imgSiftMatchLineH = '%s/img%06d_stg%d_m';         % asift match lines
%TODO vertical match
conf.imgSiftMatchHori  = '%s/img%06d_stg%d_mh';        % asift greyscal match lines horizontal alignment
conf.imgSiftMatchVert  = '%s/img%06d_stg%d_mv';        % asift greyscal match lines vertical alignment
conf.imgSiftHomoFile   = '%s/img%06d_stg%d_xh_t%1.1f'; % homography fir match lines
conf.asiftMatchRaw     = '%sasift_img%06d_stg%d';      % asift match point coordinates (text file)
conf.asiftRegionFile   = '%s/asift_match_%06d.mat';    % .mat containing structure with region and corresponding asifts points


% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.sampleSize  = length(consts.useNdx); 
conf.sampleStages = [5]; % use [5 4 3 2 1] to loop over all stages segmentations
conf.imgGap = 20; % size of gap between the images in pixels
conf.juncMarker = 'oy';
conf.siftMarker = 'oc';
conf.markerSize = 3;
conf.startFromImgID = 384;
     
for sampleStage = conf.sampleStages
for setInd = 1:length(consts.matchImgId)

idSet = consts.matchImgId{setInd};    % set of image IDs for cuurent match group (from consts.matchImgId cell array)
imgNum = length(idSet);               % number of images in chosen set

% create subdir according to Sample size and Stage chosen to keep each experiment separately
% matchDir = sprintf('%ssample%d_%d/stage_%d/', consts.matchDir, floor(idSet(1)/10^5),conf.sampleSize, sampleStage);
matchDir = sprintf('%ssample%d_000/stage_%d/', consts.matchDir, floor(idSet(1)/10^5),sampleStage); %TEMP, better use above to keep seperate dirs
if exist(matchDir, 'dir')~=7
     mkdir(matchDir);
end

if any(idSet < conf.startFromImgID); continue; end;  % FORCE TO SKIP FIRST N IMAGES (N == conf.startFromImgID);
if imgNum < 2 || imgNum > 2;         continue; end;  % CONSIDERRING ONLY PAIRS OF IMAGES, NO TRIPLETS NOR MORE (BY NOW)
% CHECK CONSISTENCY BETWEEEN SELECTED IMAGES AND USED IMAGES (const.useNdx & consts.matchImgId == naz_generate_match_ndxs)
if ~( any(consts.useNdx==idSet(1)) ||  any(consts.useNdx==idSet(2)) );  continue; end;

fprintf('\nProcessing image pairs: %d<->%d, stage %d ...\n', idSet(1), idSet(2), sampleStage); 

imRgb = cell(imgNum,1); % RGB images
im  = cell(imgNum, 1);   % greyscale images 
xj  = cell(imgNum, 1);   % junction x-coordinates
yj  = cell(imgNum, 1);   % junction y-coordinates
D   = cell(imgNum, 1);   % sift descriptors

juncsIm = cell(imgNum, 1);
edgesIm = cell(imgNum, 1);

for i = 1:imgNum 
    load(sprintf(consts.imageRgbFilename,  idSet(i)), 'imgRgb');
    if sampleStage == 0
        load(sprintf(consts.watershedFilename, idSet(i)));
    else % regardless the configuratoin, classifier trained on all images (1449) will be used
        load(sprintf(consts.boundaryInfoPostMerge, 1449, params.seg.featureSet, sampleStage, idSet(i)), 'boundaryInfo');
    end
    edgesIm{i} = boundaryInfo.edges.fragments;     % edges watershed\merging segmentation
    juncsIm{i} = boundaryInfo.junctions.position;  % round means top-right pixel from junction point
    imgRgb = double(imgRgb)/255; % im2double
    
    %MARKER rounding!
    % extracting image and its junctions coordinates
    imRgb{i} = imgRgb;
    im{i} = rgb2gray(imRgb{i});
    xj{i} = round(juncsIm{i}(:,1));
    yj{i} = round(juncsIm{i}(:,2));
    conf.sizeX = size(imRgb{i},1);
    conf.sizeY = size(imRgb{i},2);
    
    %%% DEBUG %%%
    binfo{i} = boundaryInfo; %#ok<SAGROW>
    % bi = binfo{i};
    % naz_test;    
end
clear imgRgb boundaryInfo;
conf.imgGapStub = 0*ones(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity
pairedImRgb = [imRgb{1} conf.imgGapStub imRgb{2}];

% show images
% --------------------------------------------------
filename = sprintf([conf.imgFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file') || conf.overwrite_image

    h_img = figure('Visible', 'off'); % params.debug_visible);
    imshow(pairedImRgb); axis image; axis off; title(sprintf('Images #%d #%d', idSet(1), idSet(2)));    
    if params.debug; print(h_img, '-dpng', sprintf(filename, matchDir, idSet(1), sampleStage) ); end;

end
clear filename;

% show junction points
% --------------------------------------------------
filename = sprintf([conf.imgJuncFile '.png'], matchDir, idSet(1), sampleStage);
if ~exist(filename, 'file') || conf.overwrite_image

    h_pnts = figure('Name',['Junctions, stage ' num2str(sampleStage)],'Visible', 'off');% params .debug_visible);
    imshow(pairedImRgb); axis image; axis off; hold on; 
    title(sprintf('Images #%d #%d, stage%d', idSet(1), idSet(2), sampleStage)); 
    plot(juncsIm{1}(:,1),juncsIm{1}(:,2), conf.juncMarker, 'MarkerSize', conf.markerSize);
    plot(juncsIm{2}(:,1)+conf.sizeY+conf.imgGap, juncsIm{2}(:,2), conf.juncMarker, 'MarkerSize', conf.markerSize);
    naz_plot_paired_edges(edgesIm, conf);
    
    %if params.debug;      print(h_pnts, '-dpng', sprintf('%s/img%06d_b.png', matchDir, idSet(1)) ); end;
    if params.debug;     saveas(h_pnts, filename, 'png'); end
    if params.debug_fig; saveas(h_pnts, filename, 'fig'); end

end
clear filename shift i k;

% applying a-sift to 2 images
% ------------- ------------------------------------
a.imIn{1} = 'imIn1.png';
a.imIn{2} = 'imIn2.png';
a.imgOutVert = 'imgOutVert.png';
a.imgOutHori = 'imgOutHori.png';
a.matchings = 'matchings.txt';
a.keys{1} = 'keys1.txt';
a.keys{2} = 'keys2.txt';
a.flag_resize = 0;

aiftfilepath = sprintf([conf.asiftMatchRaw '.txt'], matchDir, idSet(1), sampleStage);
if ~exist(aiftfilepath, 'file') || conf.overwrite_asift
    asiftdir = [consts.datasetDir 'indoor_scene_seg_sup/ASIFT/'];
    cd(asiftdir);
    imwrite(im{1}, a.imIn{1}, 'png');
    imwrite(im{2}, a.imIn{2}, 'png');
    demo_ASIFT(a.imIn{1}, a.imIn{2}, a.imgOutVert, a.imgOutHori, a.matchings, a.keys{1}, a.keys{2}, a.flag_resize);
    cd('..');
    % copying greyscal match images
    copyfile([asiftdir a.imgOutHori], ...
             sprintf([conf.imgSiftMatchHori '.png'], matchDir, idSet(1), sampleStage), 'f');
    copyfile([asiftdir a.imgOutVert], ...
             sprintf([conf.imgSiftMatchVert '.png'], matchDir, idSet(1), sampleStage), 'f');

    % copying match coordinates txt files    
    copyfile([asiftdir a.matchings], aiftfilepath, 'f');
else
    fprintf('Skipping asift computation, matches pre-computed in: %s\n', aiftfilepath);
end

% parse asift matching.txt 
% --------------------------
fid = fopen(aiftfilepath);
strline = regexp(fgetl(fid), ' ', 'split');
num_matches  = str2double(strline{1});
xat = zeros(1,num_matches); yat = xat; xbt = xat; ybt = xat;
for i = 1:num_matches
    strline = regexp(fgetl(fid), ' ', 'split');
    xat(i) = str2double(strline{1});
    yat(i) = str2double(strline{3});    
    xbt(i) = str2double(strline{5});
    ybt(i) = str2double(strline{7});    
end
fclose(fid);

clear aiftfilepath fid strline num_matches asiftdir;
%TODO ~exist (its legacy from previous sift algorithm)
% % filename = sprintf([conf.imgSiftHomoFile '.png'], matchDir, idSet(1), sampleStage, t, rt);
% % if ~exist(filename, 'file')


% show all tentative matches
% ----------------------------------------------
h_match = figure('Name','A-SIFT matches','Visible',params.debug_visible);
imshow(pairedImRgb); axis image; axis off; hold on;
title( sprintf('Tentative correspondences: img #%d #%d, stage%d', idSet(1), idSet(2), sampleStage) );
naz_plot_paired_edges(edgesIm, conf);
plot(xat,yat,                       conf.siftMarker, 'MarkerSize', conf.markerSize);
plot(xbt+conf.sizeY+conf.imgGap,ybt,conf.siftMarker, 'MarkerSize', conf.markerSize);

if params.debug;     saveas(h_match, sprintf([conf.imgSiftMatchPoint '.png'], matchDir, idSet(1), sampleStage), 'png'); end
if params.debug_fig; saveas(h_match, sprintf([conf.imgSiftMatchPoint '.fig'], matchDir, idSet(1), sampleStage), 'fig'); end

% plot matching line in separate figure (horizontal and [todo]vertical)
h_match_h = figure('Name','A-SIFT matches','Visible',params.debug_visible);
imshow(pairedImRgb); axis image; axis off; hold on;
title( sprintf('Tentative correspondences: img #%d #%d, stage%d', idSet(1), idSet(2), sampleStage) );
naz_plot_paired_edges(edgesIm, conf);
plot(xat,yat,                       conf.siftMarker, 'MarkerSize', conf.markerSize);
plot(xbt+conf.sizeY+conf.imgGap,ybt,conf.siftMarker, 'MarkerSize', conf.markerSize);
hl = line([xat; xbt+conf.sizeY+conf.imgGap],[yat; ybt],'color','g');

if params.debug;     saveas(h_match_h, sprintf([conf.imgSiftMatchLineH '.png'], matchDir, idSet(1), sampleStage), 'png'); end
if params.debug_fig; saveas(h_match_h, sprintf([conf.imgSiftMatchLineH '.fig'], matchDir, idSet(1), sampleStage), 'fig'); end

EXC = false;
% Robustly fit homography (removing bad matches)
% ----------------------------------------------
% Specify the inlier threshold (in noramlized image co-ordinates)
% TODO iterating through vector 'thresholds' saves only images in different files, but actual .mat data is overwritten on every iteration
thresholds = [0.5];
xah = cell(length(thresholds),1); yah = xah; xbh = xah; ybh = xah;
for i = 1:length(thresholds)
    t = thresholds(i);
    fprintf('Fitting homography with threshold = %1.1f...', t);
try
    [Hab, inliers] = ransacfithomography([xat; yat], [xbt; ybt], t);
    xah{i} = xat(inliers); yah{i} = yat(inliers); xbh{i} = xbt(inliers); ybh{i} = ybt(inliers);
    
    h_homo = figure('Name',['Fitting homography, threshold = ' num2str(t)],'Visible', params.debug_visible);
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Homography (ransac inliers): img #%d #%d, stg%d (t=%1.1f)', idSet(1), idSet(2), sampleStage, t));
    naz_plot_paired_edges(edgesIm, conf);
    plot(xah{i},yah{i},                        conf.siftMarker, 'MarkerSize', conf.markerSize);
    plot(xbh{i}+conf.sizeY+conf.imgGap,ybh{i}, conf.siftMarker, 'MarkerSize', conf.markerSize);
    %plot fitted matches
    hl = line([xah{i}; xbh{i}+conf.sizeY+conf.imgGap],[yah{i}; ybh{i}],'color','g');

    if params.debug;     saveas(h_homo, sprintf([conf.imgSiftHomoFile '.png'], matchDir, idSet(1), sampleStage,t), 'png'); end
    if params.debug_fig; saveas(h_homo, sprintf([conf.imgSiftHomoFile '.png'], matchDir, idSet(1), sampleStage,t), 'fig'); end
catch ME
    fprintf('Exception in homography: %s\n', ME.message);
    EXC = true;
    continue;
end
end

if EXC; continue; end;
%---------------------------------------------------
% group regions matches, grouping matrix Rx3, where R number of regionds
fprintf('Matching regions...\n');
ti = 1; % just pick index in 'thresholds' array, if only one threshold then put 1.

% set(0,'DefaultFigureWindowStyle','docked');
% h_reg = figure('Name','Region match','Visible',params.debug_visible); clf;
% imshow(pairedImRgb); axis image; axis off; hold on;
% title(sprintf('Region matching (based on asift): img #%d #%d', idSet(1), idSet(2)));
% naz_plot_paired_edges(edgesIm, conf);
% plot(xh{1},yh{1},                        '+', 'MarkerSize', conf.markerSize);
% plot(xh{2}+conf.sizeY+conf.imgGap,yh{2}, '+', 'MarkerSize', conf.markerSize);

for i = 1:2 % i==1 iterting through left image, i==2 - right image
if i==1
    xh{i} = xah{ti}; yh{i} = yah{ti};
else 
    xh{i} = xbh{ti}; yh{i} = ybh{ti};  
end
asiftRegions.bndrInfo{i} = binfo{i};          % boundaryInfo for left  image
asiftRegions.id(i) = idSet(i);                % image ID
asiftRegions.asiftInd{i} = [xh{i}; yh{i}];    % coordinates for matched points NB! Keep 2xN size for the plot function

% it has size of RxS, where R is number of regions, and S - number of total sift point in the image;
% matrix contains matched points per region as exhaustive logical vector
% for all the points (true - point belong to a region, false - does not belong
asiftRegions.region2ind{i} = false( length(unique(binfo{i}.imgRegions)), length(xh{i}) );  

%TODO
asiftRegions.regionContour{i} = cell(length(unique(binfo{i}.imgRegions)),1);
% regionMatch.ind2region{i}  = zeros( length(xh{i}), 1);  % TODO region corresponding to matched point

    
    
shift = 0; if i==2; shift = conf.sizeY+conf.imgGap; end;
regions = asiftRegions.bndrInfo{i}.imgRegions;
for r = unique(regions)';
    % fprintf('region = %d\n', r);
    %---------------
    % padding & removing outliers, to supress odd results from 'contour' function (bug?)
%     reg = padarray(regions==r, [1 1]);
%     cont_ini = contourc(double(reg),1);   % inital values botain from 'contour' func, with odd points <1
%     [~, cont_outlier_c] = find(cont_ini<1);
%     cont_c = logical(1:length(cont_ini));
%     cont_c(cont_outlier_c) = false;
%     cont = cont_ini(:,cont_c) - 1; % -1 to compensate coordinates shift because of padding
%     clear cont_outlier cont_outlier_c cont_c;
    %---------------
    %NEW FIX avoids akward diagonal lines (as was in presentation)
    reg = regions==r;
    %padarray(reg, [1 1]);
    [cont labels] = bwboundaries(reg, 8, 'noholes');
    lens = cellfun(@length, cont);
    [~, ind] = max(lens);
    cont = flipud(cont{ind}');

    in = inpolygon( xh{i},yh{i}, cont(1,:),cont(2,:) );
    asiftRegions.region2ind{i}(r,:) = in;
    %TODO  regionMatch.ind2region{i}(  
    annyy = any(in);
    if annyy
        hc = plot(cont(1,:)+shift, cont(2,:), 'b');
        hr = plot(xh{i}(in)+shift, yh{i}(in), 'mo', xh{i}(~in)+shift, yh{i}(~in), 'c+');
        delete(hr);
        delete(hc);
    end
end
end
save(sprintf(conf.asiftRegionFile, matchDir, idSet(1)), 'asiftRegions'); 


% pause;
close all;
end

end

