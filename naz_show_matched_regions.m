clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 
warning off all;

iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
consts.matchedDir   = [consts.datasetDir 'asift_matched/'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'on';   
params.debug_fig = false;

% Must replace below 'consts' field when using A-SIFT
conf.overwrite_image = true;
conf.imgFile = '%s/img%06d_stg%d_a';
conf.asiftMatchFile    = '%sasift_img%06d_stg%d';    % asift match point coordinates (text file)
conf.regionMatchFile   = '%s/region_img%06d.mat';     % .mat containing structure with region matching

% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.imgGap = 20; % size of gap between the images
conf.juncMarker = 'oy';
conf.siftMarker = 'oc';
conf.markerSize = 2.5;
conf.startFromImgID = 0;
conf.color = {'m', 'y', [0.06 0.9 0.4], 'w', 'k'};

matchDir = [consts.matchedDir 'sample0_1449/stage_5/'];
matlist = dir([matchDir '*.mat']);
matlist = {matlist.name};
% ===============================================================================================

for matfile = matlist
    matfile = matfile{1}; %#ok<FXSET>
    load([matchDir matfile]);
    rm = regionMatch; clear regionMatch;
    if any(rm.id < conf.startFromImgID); continue; end;  % FORCE TO SKIP FIRST N IMAGES (N == conf.startFromImgID);
    for i =1:2
        load(sprintf(consts.imageRgbFilename, rm.id(i)), 'imgRgb'); 
        imgRgb = double(imgRgb)/255; % im2double
        imRgb{i} = imgRgb; clear imgRgb; %#ok<SAGROW>
        edgesIm{i} = rm.bndrInfo{i}.edges.fragments;      %#ok<SAGROW>
        xh{i} = rm.asiftInd{i}(1,:); %#ok<SAGROW>
        yh{i} = rm.asiftInd{i}(2,:); %#ok<SAGROW>
    end

    conf.sizeX = size(imRgb{1},1); conf.sizeY = size(imRgb{1},2);    
    conf.imgGapStub = 0*ones(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
    pairedImRgb = [imRgb{1} conf.imgGapStub imRgb{2}];

    set(0,'DefaultFigureWindowStyle','docked');
    h_reg = figure('Name','Region match','Visible',params.debug_visible); clf;
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Region matching (based on asift): img #%d #%d', rm.id(1), rm.id(2)));
    naz_plot_paired_edges(edgesIm, conf, 'horizontal');
    plot(xh{1},yh{1},                        '+', 'MarkerSize', conf.markerSize);
    plot(xh{2}+conf.sizeY+conf.imgGap,yh{2}, '+', 'MarkerSize', conf.markerSize);
    

    shift = conf.sizeY+conf.imgGap;
    %NB! iterating through left image regions only
    for r = 1:size(rm.region2ind{1},1)
        in0 = rm.region2ind{1}(r,:);     
        temp = any(in0);
        if ~(any(in0)); continue; end;   % skip if region does not contain asift point
        
        % LEFT IMAGE
        regions1 = rm.bndrInfo{1}.imgRegions;
        cont = naz_region_contour(regions1==r);
        hp_0 = plot(xh{1}(in0), yh{1}(in0), 'co', 'MarkerSize', conf.markerSize);
        hc_0 = plot(cont(1,:), cont(2,:), 'g');%, 'MarkerSize', conf.markerSize);
        
        iter = 1;
        in2t = {};
        hp_1 = {}; hc_1 = {}; hp_2 = {}; hc_2 = {}; hl = {};
        for r2 = 1:size(rm.region2ind{2},1)
            in2 = rm.region2ind{2}(r2,:); 
            if any(in & in2)  % if asift points in 'r' region (from left image) is present in 'r2' region (right image)
                regions2 = rm.bndrInfo{2}.imgRegions;
                cont2 = naz_region_contour(regions2==r2);
                in_temp = inpolygon( xh{2},yh{2}, cont2(1,:),cont2(2,:) ); % re-find the inpolygon
                in2t = in2t | in_temp;
                hc_2{iter} = plot(cont2(1,:)+shift, cont2(2,:), 'g');
                hr_2 = plot(xh{2}(in2t)+shift, yh{2}(in2t), 'co', 'MarkerSize', 3);

                iter = iter + 1;
            end            
        end
%        hl1 = line([xh{1}(in2t); xh{2}(in2t)+shift],[yh{1}(in2t); yh{2}(in2t)],'color','m', 'LineWidth', 1.2);
%        hl2 = line([xh{1}(in2t); xh{2}(in2t)+shift],[yh{1}(in2t); yh{2}(in2t)],'color','m', 'LineWidth', 1.2);

         
        ht1  = text(20, 15, sprintf('%d', sum(in)), 'Color', [0 0 .8], 'FontSize', 16,  'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        ht2  = text(20+shift, 15, sprintf('%d', sum(in2t)), 'Color', [0 0 .8], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        ht22 = text(20+shift+50, 15, sprintf('%d', hc_len-1), 'Color', [0.8 0.2 0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

        
        delete(hr_1);delete(hc_1);
        delete(hr_2);        
        for ii = 1:length(hc_2); delete(hc_2{ii}); end;
        delete(hl);
        delete(ht1); delete(ht2); delete(ht22);
        clear iter hr_1 hr_2 hc_1 hc_2 hl;
        
        % manually run to see ground-truth during debuggin.
        if false
        set(0,'DefaultFigureWindowStyle','normal');
        hh = figure;
        load(sprintf(consts.instanceLabelsFilename, rm.id(1)), 'imgInstanceLabels');
        subplot(121), imshow(imgInstanceLabels,[]), axis image, axis off, colormap jet;
        load(sprintf(consts.instanceLabelsFilename, rm.id(2)), 'imgInstanceLabels');
        subplot(122), imshow(imgInstanceLabels,[]), axis image, axis off, colormap jet;
        close(hh);
        end
        
    end
    
end