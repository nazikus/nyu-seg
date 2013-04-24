%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% THIS SCRIPT RUNS ONLY IF BOUNDARY CLASSIFIER IS TRAINED ALREADY %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
addpath(genpath('.\'));
Consts;
Params;

params.debug_visible = 'on';   % doesn't work because seg2framents.m loads Params.m again
params.overwrite_feat = true;  %FIXIT is it necessary to re-compute features?
params.seg.featureSet = consts.BFT_RGBD;
OVERWRITE = true;
segtestNdxs = consts.useNdx;  % [3 40 47];

fprintf('Starting segmentation for images %s.\n', sprintf('#%d ', segtestNdxs));
load(consts.splitsPath, 'trainNdxs'); % Load the train/test split.
global segtestDir;
segtestDir = [consts.datasetDir 'naz_seg/'];
if ~exist(segtestDir, 'dir'); mkdir(segtestDir); end
fprintf('Results are saved to: %s\n', segtestDir);



for ii = segtestNdxs
    fprintf('Segmenting image #%d\n', ii);
    for stage = 1 : params.seg.numStages
        
        % loading classifiers for current stage
        fprintf('Loading classifiers for stage #%d... ', stage);
        boundaryClassifierFilename = sprintf(consts.boundaryClassifierFilename, params.seg.featureSet, stage);
        boundaryClassifierFilename = sprintf('%s_s%d.mat', boundaryClassifierFilename(1:end-4), length(consts.useNdx) );

        assert(exist(boundaryClassifierFilename, 'file')~=0, '[Assert error] file does not exist: %s\n', boundaryClassifierFilename);
        load(boundaryClassifierFilename, 'classifier'); 
        fprintf('done.\n');

        if consts.useImages(ii) && any(trainNdxs==ii)
          fprintf('\nImage #%4d was used for training classifiers!\n\n', ii);
        end
        
        % Extracting boundary classifier features (regardles, because 'params.overwrite_feat = true')
        extract_boundary_classifier_features_and_labels(stage, params, ii, ii);    

        % Loading merging data
        fprintf('Merging regions (Image %d/%d, stage %d)... ', ii, consts.numImages, stage);
        outFilename = sprintf(consts.boundaryInfoPostMerge, params.seg.featureSet, stage, ii);
        if exist(outFilename, 'file') && ~OVERWRITE
            fprintf(' already exists, overwrite=false.\n');
            continue;
        end

        %FIXIT check if watershed is precomputed, if not - compute it
        if stage == 1
          boundaryInfoFilename = sprintf(consts.watershedFilename, ii);
        else
          boundaryInfoFilename = sprintf(consts.boundaryInfoPostMerge, params.seg.featureSet, stage-1, ii);
        end
        load(boundaryInfoFilename, 'boundaryInfo');
        load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
        load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
        load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
        load(sprintf(consts.boundaryFeaturesFilename, params.seg.featureSet, stage, ii), 'boundaryFeatures');

        [~, instanceLabels] = get_labels_from_instances(boundaryInfo.imgRegions, imgObjectLabels, imgInstanceLabels);

        % Merging
        result = merge_regions(boundaryInfo, boundaryFeatures, classifier, stage, params);
        boundaryInfo = update_boundary_info(boundaryInfo, result, imgRgb, ii, stage);
        save(outFilename, 'boundaryInfo');

        fprintf(' Done.\n');
    end  
    fprintf('\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n');
    fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n');
end

clearvars -global segtestDir;
