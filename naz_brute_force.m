clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 
warning off all
iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
consts.matchedDir   = [consts.datasetDir 'asift_matched/'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'off';   

conf.startFromImgID = 1;%173;
conf.imgGap = 20; % size of gap between the images
conf.juncMarker = 'oy';
conf.siftMarker = 'oc';
conf.markerSize = 2;
conf.lineWdith = 0.75;
conf.color = {[0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4]};
desat = @(rgb) hsv2rgb(rgb2hsv(rgb)*[1 0 0; 0 0.6 0; 0 0 1]); % de-saturate - decrease RGB saturation
conf.method = 'bruteforce';

% Must replace below 'consts' field when using A-SIFT
conf.overwrite_image = true;
conf.penFile = '%s/img%06d';

%conf.asiftRegionFile   = '%s/asift_match_%06d.mat';  % .mat containing structure with region matching
matchDir = [consts.matchedDir 'sample0_1449/stage_5/'];
matlist = dir([matchDir '*.mat']);
matlist = {matlist.name};
% ===============================================================================================
sumUnweightedScoreIni = 0;
sumWeightedScoreIni = 0;
sumUnweightedScorePen = 0;
sumWeightedScorePen = 0;
count = 0;
for matfile = matlist
matfile = matfile{1}; %#ok<FXSET>
load([matchDir matfile]);
rm = regionMatch; clear regionMatch;  %ar = asiftRegions; clear asiftRegions;
matchedFrags{1} = matchedFragsLR;
matchedFrags{2} = matchedFragsRL;
if any(rm.id < conf.startFromImgID); continue; end;  % FORCE TO SKIP FIRST N IMAGES (N - conf.startFromImgID);
for i =1:2
    load(sprintf(consts.imageRgbFilename, rm.id(i)), 'imgRgb'); 
    imgRgb = double(imgRgb)/255; % im2double
    imRgb{i} = imgRgb; clear imgRgb; %#ok<SAGROW>
%     edgesIm{i} = rm.bndrInfo{i}.edges.fragments;      %#ok<SAGROW>
%     xa{i} = rm.asiftInd{i}(1,:); %#ok<SAGROW> % A-SIFT point-x %#ok<SAGROW>
%     ya{i} = rm.asiftInd{i}(2,:); %#ok<SAGROW> % A-SIFT point-y %#ok<SAGROW>
end

% LR=1 - Left-to-Right matching orientation, LR=2 - Right-to-Left orientation
for LR = 1:2
    I = LR;
    %I2 = mod(LR,2)+1;  % I2 == 2 if LR == 1, and I2 == 1 if LR == 2
    conf.sizeX = size(imRgb{I},1); conf.sizeY = size(imRgb{I},2);    
    conf.imgGapStub = zeros(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
    pairedImRgb = [imRgb{I} conf.imgGapStub imRgb{I}];
    shift = conf.sizeY + conf.imgGap;

    penBnd = naz_remove_fragments(rm.bndrInfo{I}, matchedFrags{I});
    edgesIm{1} = rm.bndrInfo{I}.edges.fragments;
    edgesIm{2} = penBnd.edges.fragments;
    
    fprintf('Processing image # %d\n', rm.id(I));
    load(sprintf(consts.objectLabelsFilename, rm.id(I)), 'imgObjectLabels');
    load(sprintf(consts.instanceLabelsFilename, rm.id(I)), 'imgInstanceLabels');
    imgRegionsTrue = get_regions_from_labels(imgObjectLabels, imgInstanceLabels); clear imgObjectLabels imgInstanceLabels;
    % ws - weighted score, us - unweighted score, Ini - initial segmentation, Pen - penelaized segmentation
    [wsIni, usIni] = evaluate_segmentation(imgRegionsTrue, rm.bndrInfo{I}.imgRegions);
    [wsPen, usPen] = evaluate_segmentation(imgRegionsTrue, penBnd.imgRegions);
    sumWeightedScoreIni = sumWeightedScoreIni + wsIni;
    sumUnweightedScoreIni = sumUnweightedScoreIni + usIni;
    sumWeightedScorePen = sumWeightedScorePen + wsPen;
    sumUnweightedScorePen = sumUnweightedScorePen + usPen;
    
    set(0,'DefaultFigureWindowStyle','normal');
    h_fig = figure('Name','Region match','Visible',params.debug_visible); clf;
    subplot(211),
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Penalization (%s): img %d\n%.2f/%.2f vs %.2f/%.2f', ...
       conf.method, rm.id(I), wsIni, usIni, wsPen, usPen));
    naz_plot_paired_edges(edgesIm, conf);
    
    for k=1:size(matchedFrags{I},1)
    try
        naz_plot_fragments(matchedFrags{I}{k,2}, rm.bndrInfo{I}, 0, conf.color{k}, conf.lineWdith );
    catch ME
       fprintf('#### ERROR: %s; (Iteration %d)\n', ME.message, k);
    end
    end
    
    subplot(212),
    image([0, conf.sizeY],[0, conf.sizeX], imgRegionsTrue, 'CDataMapping','scaled');
    title('ground truth'), axis image, axis off;
    % image([conf.sizeY+conf.imgGap, conf.sizeY*2+conf.imgGap],[0, conf.sizeX], imgRegionsTrue, 'CDataMapping','scaled');
    saveas(h_fig, sprintf([conf.penFile '.png'], [matchDir conf.method], rm.id(I)), 'png');
    close(h_fig);
    count = count + 1;
    fprintf('w ini av = %.2f, uw ini av = %.2f\nw pen av = %.2f, uw pen av = %.2f\n\n', ...
        sumWeightedScoreIni/count, sumUnweightedScoreIni/count, sumWeightedScorePen/count, sumUnweightedScorePen/count);
end
    
end
