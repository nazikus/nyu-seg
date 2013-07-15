clear; clc; close all;  addpath(genpath('.\')); % ind_ = @(A,r,c) A(r,c); 
warning off all
profile clear
profile off
%profile -memory on
diary on;

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
conf.asiftRegionFile   = '%s/asift_match_%06d.mat';  % .mat containing structure with region matching
conf.fragPenalizeTotal = '%s/img%06d_xp';            % image pair depicting which fragment to be penalized and the ground truth
conf.fragPenalize      = '%s/img%06d_%02d_xp';       % image pair depicting which fragment to be penalized and the ground truth, per iteration
conf.regionLLL         = '%s/img%06d_%02d_xl';       % actual region matching
conf.regionMMM         = '%s/img%06d_%02d_xm';       % actual region matching

% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
conf.startFromImgID = 0;
conf.imgGap = 20; % size of gap between the images
conf.juncMarker = 'oy';
conf.siftMarker = 'oc';
conf.markerSize = 2;
conf.lineWdith = 0.75;
conf.color = {[0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4],    [1 0 1], [1 1 0], [0.06 0.9 0.4], [0.65 0.7 0.43], [0 1 0], [1 1 0], [0 1 1], [1 0 1],[0.65 0.7 0.43], [0.06 0.9 0.4]};
desat = @(rgb) hsv2rgb(rgb2hsv(rgb)*[1 0 0; 0 0.6 0; 0 0 1]); % de-saturate - decrease RGB saturation

matchDir = [consts.matchedDir 'sample0_001/stage_5/'];
matlist = dir([matchDir '*.mat']);
matlist = {matlist.name};
% ===============================================================================================

for matfile = matlist
    matfile = matfile{1}; %#ok<FXSET>
    load([matchDir matfile]);
    % rm = regionMatch; clear regionMatch;  
    rm = asiftRegions; clear asiftRegions;
    
    if any(rm.id < conf.startFromImgID); continue; end;  % FORCE TO SKIP FIRST N IMAGES (N - conf.startFromImgID);
    for i =1:2
        load(sprintf(consts.imageRgbFilename, rm.id(i)), 'imgRgb'); 
        imgRgb = double(imgRgb)/255; % im2double
        imRgb{i} = imgRgb; clear imgRgb; %#ok<SAGROW>
        edgesIm{i} = rm.bndrInfo{i}.edges.fragments;      %#ok<SAGROW>
        xa{i} = rm.asiftInd{i}(1,:); %#ok<SAGROW> % A-SIFT point-x %#ok<SAGROW>
        ya{i} = rm.asiftInd{i}(2,:); %#ok<SAGROW> % A-SIFT point-y %#ok<SAGROW>
    end  
    
    % LR=1 - Left-to-Right matching, LR=2 - Right-to-Left matching 
    for LR = 1:2
    tic;
    I1 = LR;
    I2 = mod(LR,2)+1;  % I2 == 2 if LR == 1, and I2 == 1 if LR == 2
    fprintf('Processing image pair (%d, %d), N(asift)==%d\n', rm.id(I1), rm.id(I2), length(rm.asiftInd{LR}));
    
    if length(rm.asiftInd{1})<400; % skip images with more then N asift matches
        fprintf('Skipping...\n\n');
        continue; 
    end; 

    %if skip already computed
    
    edgesIm = {edgesIm{I1}; edgesIm{I2}};
    conf.sizeX = size(imRgb{I1},1); conf.sizeY = size(imRgb{I2},2);    
    conf.imgGapStub = zeros(conf.sizeX, conf.imgGap, 3); % 1 == maximum intensity (255)
    pairedImRgb = [imRgb{I1} conf.imgGapStub imRgb{I2}];
    shift = conf.sizeY + conf.imgGap;
    
    set(0,'DefaultFigureWindowStyle','docked');
    h_fig = figure('Name','Region match','Visible',params.debug_visible); clf;
    imshow(pairedImRgb); axis image; axis off; hold on;
    title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(I1), rm.id(I2)));
    naz_plot_paired_edges(edgesIm, conf);

    %NB! 
    % if LR ==1
    % iterating through right image regions, to penalize edges in the left image
    % suffix R deontes right image, L - left image
    % in -  bool values, which asift points are within the currently iterated region
    % h - handler, c - contour, a - asift, l - line
    % if LR ==2, then left-right is opposite!
    load(sprintf(consts.objectLabelsFilename, rm.id(I1)), 'imgObjectLabels');
    load(sprintf(consts.instanceLabelsFilename, rm.id(I1)), 'imgInstanceLabels');
    imgRegionsTrue = get_regions_from_labels(imgObjectLabels, imgInstanceLabels); clear imgObjectLabels imgInstanceLabels;

    regCentriods = regionprops(rm.bndrInfo{I1}.imgRegions, 'Centroid');
    matchedRegions = {}; matchedFrags = {};
    it = 1; % iterator over regions on the right image (that have at least 1 match with left image)
    for rR = 1:size(rm.region2ind{I2},1)
        inR = rm.region2ind{I2}(rR,:);     
        if ~(any(inR)); continue; end;   % skip if region does not contain asift point
        haL_all = plot(xa{I1},ya{I1},'ob', 'MarkerSize', conf.markerSize);
        haR_all = plot(xa{I2}+conf.sizeY+conf.imgGap,ya{I2}, 'ob', 'MarkerSize', conf.markerSize);
        
        regionsR = rm.bndrInfo{I2}.imgRegions;
        contR = naz_region_contour(regionsR==rR);
        hc_R = plot(contR(1,:)+shift, contR(2,:), 'Color', conf.color{1}, 'LineWidth', conf.lineWdith );
        ha_R = plot(xa{I2}(inR)+shift, ya{I2}(inR), 'o', 'Color', conf.color{3}, 'MarkerSize', conf.markerSize);
        
        iter = 1; % iterator over regions from left image matching to a single regions on the right image
        inLt={}; % total in left image, but each cell belongs to seperate region
        hc_L={}; ha_L={}; hl={};
        matchedRegions{it} = []; %#ok<SAGROW>
        for rL = 1:size(rm.region2ind{I1},1)
            inL = rm.region2ind{I1}(rL,:); 
            if any(inL & inR)  % if asift points in rL region (from left image) is present in rR region (right image)
                regionsL = rm.bndrInfo{I1}.imgRegions;
                contL = naz_region_contour(regionsL==rL);
                inLt{iter} = inL; %#ok<SAGROW>
                hc_L{iter} = plot(contL(1,:), contL(2,:), 'Color', conf.color{ iter }, 'LineWidth', conf.lineWdith ); %#ok<SAGROW>
                ha_L{iter} = plot(xa{I1}(inL&inR), ya{I1}(inL&inR), 'o', 'MarkerSize', conf.markerSize, 'Color', (conf.color{iter}) ); %#ok<SAGROW>
                hl{iter}   = line([xa{I1}(inL&inR); xa{I2}(inL&inR)+shift],[ya{I1}(inL&inR); ya{I2}(inL&inR)],'color', desat(conf.color{iter}), 'LineWidth', conf.lineWdith ); %#ok<SAGROW>
                iter = iter + 1;
                matchedRegions{it} = [matchedRegions{it} rL]; %#ok<SAGROW>
            end            
        end
        iter = iter-1;
        if iter>1
            [matchedFrags{it,1} matchedFrags{it,2}] = naz_get_fragments(matchedRegions{it}, rm.bndrInfo{I1}); %#ok<SAGROW>
            matchedFrags{it,3} = rR;  % persist which region on the right was matched
        elseif iter == 0
            continue;
        else
            matchedRegions = matchedRegions(1:end-1);
            delete(hc_L{1}, ha_L{1}, hl{1}, ha_R, hc_R, haL_all, haR_all); clear hc_L ha_L hl ha_R hc_R haL_all haR_all;
            continue;
        end
        
        fprintf('   regionR: %d (iter %d), matchesL: %d, inR: %d, inLt: %d\n', rR, it, iter, sum(inR), sum(cell2mat(inLt)));
        msg_htR = sprintf('%d', sum(inR));
        msg_htL = sprintf('%d', sum(cell2mat(inLt)));
        msg_htLr = sprintf('%d', iter);
        htL  = text(20 + 50,  20, msg_htL,  'Color', [0 0 0], 'FontSize', 13, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        htR  = text(20+shift, 20, msg_htR,  'Color', conf.color{3}, 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        htLr = text(20,       20, msg_htLr, 'Color', [1 .4  0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %TODO create dir 'penalized' if absent, so far must be created manually
        saveas(h_fig, sprintf([conf.regionLLL '.png'], [matchDir 'penalized'], rm.id(I1), it), 'png');

        for ii = 1:iter; delete(hl{ii});end; clear hl;
        saveas(h_fig, sprintf([conf.regionMMM '.png'], [matchDir 'penalized'], rm.id(I1), it), 'png');
        
        delete(ha_R, hc_R, haL_all, haR_all); clear ha_R hc_R haL_all haR_all;
        delete(htL, htR, htLr); clear htL htR htLr;
        for ii = 1:iter; delete(hc_L{ii}, ha_L{ii}); end;  clear hc_L ha_L;

        title(sprintf('Region penalization for img #%d, regionR #%d', rm.id(I1), rR));
        naz_plot_fragments(matchedFrags{it,2}, rm.bndrInfo{I1}, 0, 'c', conf.lineWdith );
        image([conf.sizeY+conf.imgGap, conf.sizeY*2+conf.imgGap],[0, conf.sizeX], imgRegionsTrue, 'CDataMapping','scaled');
        for k = 1:iter
            cd = regCentriods(matchedRegions{it}(k)).Centroid;
            text(cd(1), cd(2), num2str(sum(inLt{k})),  'Color', 'g', 'FontSize', 9, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            text(cd(1)+20, cd(2), num2str(sum(inLt{k}&inR)),  'Color', [0.9 0.9 0.9], 'FontSize', 11, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        text(20, 20, num2str(iter),  'Color', [1 .4  0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        text(20+shift, 20, num2str(sum(inR)),  'Color', [0.9 0.9 0.9], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
        saveas(h_fig, sprintf([conf.fragPenalize '.png'], [matchDir 'penalized'], rm.id(I1), it), 'png');
        it = it + 1; 
        
%         htT = text(20,       20, sprintf('%d',iter), 'Color', [1 .4  0], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         delete(htT); clear htT;
        
        % FINISHING ITERATION
        % getting back the original figure with paired images
        % delete(h_fig);
        % h_fig = imshow(pairedImRgb); axis image; axis off; hold on;
        image(pairedImRgb);
        title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(I1), rm.id(I2)));
        naz_plot_paired_edges(edgesIm, conf);
        
        % =================================================
        % manually run to see ground-truth during debuggin.
        if false
        set(0,'DefaultFigureWindowStyle','normal'); %#ok<UNRCH>
        hh = figure;
        load(sprintf(consts.objectLabelsFilename, rm.id(I1)), 'imgObjectLabels');
        subplot(121), imshow(imgObjectLabels,[]), axis image, axis off, colormap jet;
        load(sprintf(consts.objectLabelsFilename, rm.id(I2)), 'imgObjectLabels');
        subplot(122), imshow(imgObjectLabels,[]), axis image, axis off, colormap jet;
        close(hh);
        end
        
    end
    
    if ~isempty(matchedFrags)
        %%% saving (appending) fragments matched (potentially to be penalized)
        matchedFragsLR = matchedFrags; matchedFragsRL = matchedFrags;
        if LR==1
            save([matchDir matfile], 'matchedFragsLR', '-append');
        else
            save([matchDir matfile], 'matchedFragsRL', '-append');
        end
        clear matchedFragsLR matchedFragsRL;
        
        %%% plot all fragments at once in file (thos to be penalized)
        title(sprintf('Region penalization for img #%d', rm.id(I1)));
        for k=1:size(matchedFrags,1)
        try
            naz_plot_fragments(matchedFrags{k,2}, rm.bndrInfo{I1}, 0, conf.color{k}, conf.lineWdith );
        catch ME
           fprintf('#### ERROR: %s; (Iteration %d)\n', ME.message, k);
        end
        end
        image([conf.sizeY+conf.imgGap, conf.sizeY*2+conf.imgGap],[0, conf.sizeX], imgRegionsTrue, 'CDataMapping','scaled');
        % text = total number of fragments detected
        text(20+shift, 20, num2str(sum(cellfun(@(x) length(x), matchedFrags(:,1)))),  'Color', [0.9 0.9 0.9], 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        saveas(h_fig, sprintf([conf.fragPenalizeTotal '.png'], [matchDir 'penalized'], rm.id(I1)), 'png');
%        delete(h_fig);

%        h_fig = imshow(pairedImRgb); axis image; axis off; hold on;
        image(pairedImRgb);
        title(sprintf('Region matching (based on asift): img pair (%d, %d)', rm.id(I1), rm.id(I2)));
        naz_plot_paired_edges(edgesIm, conf);
        
        % apend found fragment matches to initial mat file
    end
    fprintf('\n ------------------------------------------\n');
    fprintf('Time passed: %ds\n', round(toc));
    %profile viewer;
    diary([matchDir 'penalized/' 'region_match.log']);
    
    end
end
diary off;

% back matching of asifts
% indif = (inL-inR)>0;
% if ~isempty(indif)
% hl{iter}   = line([xa{1}(indif); xa{2}(indif)+shift],[ya{1}(indif); ya{2}(indif)],'color', desat(conf.color{iter*2}), 'LineWidth', 1.2);
% hx{iter} = plot([xa{1}(indif) xa{2}(indif)+shift], [ya{1}(indif) ya{2}(indif)], 'o', 'MarkerSize', 3, 'Color', conf.color{iter*2});
% end
