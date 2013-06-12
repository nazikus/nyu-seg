clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 
warning off all;

iptsetpref('ImshowBorder','tight'); %??

Consts; Params;
consts.matchedDir   = [consts.datasetDir 'asift_matched/'];

params.seg.featureSet = consts.BFT_RGBD;
params.debug = true;
params.debug_visible = 'off';   
params.debug_fig = false;

% Must replace below 'consts' field when using A-SIFT
conf.overwrite_image = true;
conf.imgFile = '%s/img%06d_stg%d_a';
conf.asiftMatchFile    = '%sasift_img%06d_stg%d';    % asift match point coordinates (text file)
conf.regionMatchFile   = '%s/region_img%06d.mat';    % .mat containing structure with region matching
conf.fragPenalizeTotal = '%s/img%06d_xp';      % image pair depicting which fragment to be penalized and the ground truth
conf.fragPenalize      = '%s/img%06d_%02d_xp';  % image pair depicting which fragment to be penalized and the ground truth, per iteration
conf.regionLLL         = '%s/img%06d_%02d_xl';  % actual region matching
conf.regionMMM         = '%s/img%06d_%02d_xm';  % actual region matching

% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.imgGap = 20; % size of gap between the images
conf.juncMarker = 'oy';
conf.siftMarker = 'oc';
conf.markerSize = 2;
conf.lineWdith = 0.5;
conf.startFromImgID = 0;
conf.color = {[0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4]};
desat = @(rgb) hsv2rgb(rgb2hsv(rgb)*[1 0 0; 0 0.6 0; 0 0 1]); % de-saturate - decrease RGB saturation

matchDir = [consts.matchedDir 'sample0_1449/stage_5/'];
matlist = dir([matchDir '*.mat']);
matlist = {matlist.name};
% ===============================================================================================

for matfile = matlist
    matfile = matfile{1}; %#ok<FXSET>
    load([matchDir matfile]);
    rm = regionMatch; clear regionMatch;  %ar = asiftRegions; clear asiftRegions;
    if any(rm.id < conf.startFromImgID); continue; end;  % FORCE TO SKIP FIRST N IMAGES (N - conf.startFromImgID);
    for i =1:2
        load(sprintf(consts.imageRgbFilename, rm.id(i)), 'imgRgb'); 
        imgRgb = double(imgRgb)/255; % im2double
        imRgb{i} = imgRgb; clear imgRgb; %#ok<SAGROW>
        edgesIm{i} = rm.bndrInfo{i}.edges.fragments;      %#ok<SAGROW>
        xa{i} = rm.asiftInd{i}(1,:); % A-SIFT point-x %#ok<SAGROW>
        ya{i} = rm.asiftInd{i}(2,:); % A-SIFT point-y %#ok<SAGROW>
    end

    conf.sizeX = size(imRgb{1},1); conf.sizeY = size(imRgb{1},2);    
    conf.imgGapStub = 0*ones(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
    pairedImRgb = [imRgb{1} conf.imgGapStub imRgb{2}];
    shift = conf.sizeY + conf.imgGap;

    set(0,'DefaultFigureWindowStyle','docked');
    h_fig = figure('Name','Region match','Visible',params.debug_visible); clf;
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(1), rm.id(2)));
    naz_plot_paired_edges(edgesIm, conf);

    %NB! iterating through right image regions, to penalize edges in the left image
    % suffix R deontes right image, L - left image
    % in -  bool values, which points are within the currently iterated region
    % h - handler, c - contour, a - asift, l - line
    load(sprintf(consts.instanceLabelsFilename, rm.id(1)), 'imgInstanceLabels');
    regCentriods = regionprops(rm.bndrInfo{1}.imgRegions, 'Centroid');
    matchedRegions = {};
    it = 1;
    for rR = 1:size(rm.region2ind{2},1)
        inR = rm.region2ind{2}(rR,:);     
        if ~(any(inR)); continue; end;   % skip if region does not contain asift point
        haL_all = plot(xa{1},ya{1},'.b', 'MarkerSize', conf.markerSize);
        haR_all = plot(xa{2}+conf.sizeY+conf.imgGap,ya{2}, '.b', 'MarkerSize', conf.markerSize);
        
        regionsR = rm.bndrInfo{2}.imgRegions;
        contR = naz_region_contour(regionsR==rR);
        hc_R = plot(contR(1,:)+shift, contR(2,:), 'Color', conf.color{1}, 'LineWidth', conf.lineWdith );
        ha_R = plot(xa{2}(inR)+shift, ya{2}(inR), 'o', 'Color', conf.color{3}, 'MarkerSize', conf.markerSize);
        
        iter = 1; 
        inLt={}; % total in left image, but each cell belongs to seperate region
        hc_L={}; ha_L={}; hl={};
        matchedRegions{it} = []; %#ok<SAGROW>
        for rL = 1:size(rm.region2ind{1},1)
            inL = rm.region2ind{1}(rL,:); 
            if any(inL & inR)  % if asift points in rL region (from left image) is present in rR region (right image)
                regionsL = rm.bndrInfo{1}.imgRegions;
                contL = naz_region_contour(regionsL==rL);
                inLt{iter} = inL; %#ok<SAGROW>
                hc_L{iter} = plot(contL(1,:), contL(2,:), 'Color', conf.color{ iter }, 'LineWidth', conf.lineWdith ); %#ok<SAGROW>
                ha_L{iter} = plot(xa{1}(inL&inR), ya{1}(inL&inR), 'o', 'MarkerSize', conf.markerSize, 'Color', (conf.color{iter}) ); %#ok<SAGROW>
                hl{iter}   = line([xa{1}(inL&inR); xa{2}(inL&inR)+shift],[ya{1}(inL&inR); ya{2}(inL&inR)],'color', desat(conf.color{iter}), 'LineWidth', conf.lineWdith ); %#ok<SAGROW>
                iter = iter + 1;
                matchedRegions{it} = [matchedRegions{it} rL]; %#ok<SAGROW>
            end            
        end
        iter = iter-1;
        if iter>1
            [matchedFrags{it,1} matchedFrags{it,2}] = naz_get_fragments(matchedRegions{it}, rm.bndrInfo{1}); %#ok<SAGROW>
        else
            matchedRegions = matchedRegions(1:end-1);
            delete(hc_L{1}, ha_L{1}, hl{1}, ha_R, hc_R, haL_all, haR_all);
            continue;
        end
        
        fprintf('Img pair (%d, %d), regionR: %d, matchesL: %d\n', rm.id(1), rm.id(2), rR, iter);
        msg_htR = sprintf('%d', sum(inR));
        msg_htL = sprintf('%d', sum(cell2mat(inLt)));
        msg_htLr = sprintf('%d', iter);
        htL  = text(20 + 50,  20, msg_htL,  'Color', [0 0 0], 'FontSize', 13, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        htR  = text(20+shift, 20, msg_htR,  'Color', conf.color{3}, 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        htLr = text(20,       20, msg_htLr, 'Color', [1 .4  0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        saveas(h_fig, sprintf([conf.regionLLL '.png'], [matchDir 'penalized'], rm.id(1), it), 'png');

        for ii = 1:iter; delete(hl{ii});end;
        saveas(h_fig, sprintf([conf.regionMMM '.png'], [matchDir 'penalized'], rm.id(1), it), 'png');
        
        delete(ha_R, hc_R, haL_all, haR_all);
        delete(htL, htR, htLr); clear htL htR htLr;

        for ii = 1:iter; delete(hc_L{ii}, ha_L{ii}); end;  % delete(hl{ii});
        clear ha_L ha_R hc_L hc_R hl haL_all haR_all;

        title(sprintf('Region penalization for img #%d, regionR #%d', rm.id(1), rR));
        naz_plot_fragments(matchedFrags{it,2}, rm.bndrInfo{1}, 0, 'c', conf.lineWdith );
        image([conf.sizeY+conf.imgGap, conf.sizeY*2+conf.imgGap],[0, conf.sizeX], imgInstanceLabels, 'CDataMapping','scaled');
        for k = 1:iter
            cd = regCentriods(matchedRegions{it}(k)).Centroid;
            text(cd(1), cd(2), num2str(sum(inLt{k})),  'Color', 'g', 'FontSize', 9, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            text(cd(1)+20, cd(2), num2str(sum(inLt{k}&inR)),  'Color', [0.9 0.9 0.9], 'FontSize', 11, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        text(20, 20, num2str(iter),  'Color', [0.9 0.9 0.9], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        text(20+shift, 20, num2str(sum(inR)),  'Color', [0.9 0.9 0.9], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
        saveas(h_fig, sprintf([conf.fragPenalize '.png'], [matchDir 'penalized'], rm.id(1), it), 'png');
        it = it + 1; 
        
%         htT = text(20,       20, sprintf('%d',iter), 'Color', [1 .4  0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         delete(htT); clear htT;
        
        % FINISHING ITERATION
        % getting back original figure with paired images
        imshow(pairedImRgb); axis image; axis off; hold on;
        title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(1), rm.id(2)));
        naz_plot_paired_edges(edgesIm, conf);
        
        % =================================================
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
    
    %%% plot all fragments to be penalized at once in one file
        title(sprintf('Region penalization for img #%d', rm.id(1)));
    for k=1:size(matchedFrags,1)
        naz_plot_fragments(matchedFrags{k,2}, rm.bndrInfo{1}, 0, conf.color{k}, conf.lineWdith );
    end
    image([conf.sizeY+conf.imgGap, conf.sizeY*2+conf.imgGap],[0, conf.sizeX], imgInstanceLabels, 'CDataMapping','scaled');
    saveas(h_fig, sprintf([conf.fragPenalizeTotal '.png'], [matchDir 'penalized'], rm.id(1)), 'png');
    
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(1), rm.id(2)));
    naz_plot_paired_edges(edgesIm, conf);

end

% back matching of asifts
% indif = (inL-inR)>0;
% if ~isempty(indif)
% hl{iter}   = line([xa{1}(indif); xa{2}(indif)+shift],[ya{1}(indif); ya{2}(indif)],'color', desat(conf.color{iter*2}), 'LineWidth', 1.2);
% hx{iter} = plot([xa{1}(indif) xa{2}(indif)+shift], [ya{1}(indif) ya{2}(indif)], 'o', 'MarkerSize', 3, 'Color', conf.color{iter*2});
% end
