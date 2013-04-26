%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% THIS SCRIPT RUNS ONLY IF BOUNDARY CLASSIFIER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% FOR CORRESPONDING SAMPLE SIZE HAS BEEN TRAINED ALREADY %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; addpath(genpath('.\'));
Consts; Params;

OVERWRITE_E = true;
OVERWRITE_M = true;
params.seg.featureSet = consts.BFT_RGBD;
params.debug_visible = 'off';   % doesn't work because seg2framents.m loads Params.m again
sampleSize  = length(consts.useNdx); %%% NB! Do not change this line, change sample size only by changing the range of consts.useNdx!
global segtestDir; segtestDir = [consts.datasetDir consts.segmentDir consts.sampleDir];
segtestNdxs = [471]; % consts.useNdx;

fprintf('Starting segmentation for %d images: %s\n', length(segtestNdxs) ,sprintf('#%d ', segtestNdxs));
fprintf('Classifier trained on sample size: %d\n', sampleSize);

% check if boundary classifier has been already trained for given sample size
assert(exist(consts.boundaryFeaturesDir, 'dir')==7, ...
    'There is no trained classifier for %d images.\n "%s" dir does not exist\n\n', sampleSize, consts.boundaryFeaturesDir);

% create subdir for given sample size
if ~exist(segtestDir, 'dir'); mkdir(segtestDir); end

load(consts.splitsPath, 'trainNdxs'); % Load the train/test split.

for ii = segtestNdxs
    if any(trainNdxs==ii) && ii<consts.useNdx
      setname = 'Train';
    else
      setname = 'Test';
    end
    fprintf('\nSegmenting image #%d (its in %s set).', ii);
    for stage = 1 : params.seg.numStages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loading classifiers for current stage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nLoading classifiers for stage %d... ', stage);
        boundaryClassifierFilename = sprintf(consts.boundaryClassifierFilename, sampleSize, params.seg.featureSet, stage);

        assert(exist(boundaryClassifierFilename, 'file')==2, '**Assert error** file does not exist: %s\n', boundaryClassifierFilename);
        load(boundaryClassifierFilename, 'classifier'); 
        fprintf('done.\n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extracting boundary classifier features for current image in current stage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Extracting boundary features for image #%d (Stage %d)... ', ii, stage);
        boundaryFeaturesFilename  = sprintf(consts.boundaryFeaturesFilename, sampleSize, params.seg.featureSet, stage, ii);
        boundaryFeaturesFilenameSeg = strrep(boundaryFeaturesFilename, consts.boundaryDir, consts.segmentDir);
        
        if exist(boundaryFeaturesFilenameSeg, 'file') && ~OVERWRITE_E
            fprintf('skipping (exists), overwrite=false.\n'); 
            
        elseif exist(boundaryFeaturesFilename, 'file') && ~OVERWRITE_E
            fprintf('exists, copying from "%s".\n', consts.boundaryFeaturesDir);
            [status message] = copyfile(boundaryFeaturesFilename, boundaryFeaturesFilenameSeg, 'f');
            %fprintf('**DEBUG** Copy status: %d; Message: "%s" **\n', status, message);
            clear status message;
        
        else % actual extraction (similar to extract_boundary_classifier_features_and_labels.m)
            clear boundaryFeaturesFilename;
            if stage == 1
              prevBoundaryDataFilenameSeg = sprintf(consts.watershedFilename, ii);
            else
              prevBoundaryDataFilenameSeg = sprintf(consts.boundaryInfoPostMerge, length(consts.useNdx), params.seg.featureSet, stage-1, ii);
              prevBoundaryDataFilenameSeg = strrep(prevBoundaryDataFilenameSeg, consts.boundaryDir, consts.segmentDir);
            end
            load(prevBoundaryDataFilenameSeg, 'boundaryInfo');
            load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
            load(sprintf(consts.planeDataFilename, ii), 'planeData'); 
            load(sprintf(consts.watershedFilename, ii), 'pbAll');
            load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
            load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');

            %get list of each instance label for each region
            [~, instanceLabels] = get_labels_from_instances(boundaryInfo.imgRegions, imgObjectLabels, imgInstanceLabels);    
            [boundaryFeatures, boundaryLabels] = get_boundary_classifier_features(ii, imgRgb, planeData, boundaryInfo, pbAll, instanceLabels, params);
            save(boundaryFeaturesFilenameSeg, 'boundaryFeatures', 'boundaryLabels');
            fprintf('done.\n');
            clear boundaryInfo imgRgb planeData pbAll imgObjectLabels imgInstanceLabels;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loading merging data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Merging regions for image #%d (stage %d)... ', ii, stage);
        postMergeFilenameSeg = sprintf(consts.boundaryInfoPostMerge, sampleSize, params.seg.featureSet, stage, ii);
        postMergeFilenameSeg = strrep(postMergeFilenameSeg, consts.boundaryDir, consts.segmentDir);
        if exist(postMergeFilenameSeg, 'file') && ~OVERWRITE_M
            fprintf(' already exists, overwrite=false.\n');
        else

            %TODO Compute all watersheds!
            if stage == 1
              boundaryInfoFilenameSeg = sprintf(consts.watershedFilename, ii);
              source      = sprintf('%s%06d_d_fragments.png',consts.watershedDir,ii);
              destination = sprintf('%ss%d_%06d_fragments_stage0.png', segtestDir, length(consts.useNdx), ii);
              copyfile(source, destination); clear source destination;
            else
              boundaryInfoFilenameSeg = sprintf(consts.boundaryInfoPostMerge, sampleSize, params.seg.featureSet, stage-1, ii);
              boundaryInfoFilenameSeg = strrep(boundaryInfoFilenameSeg, consts.boundaryDir, consts.segmentDir);
            end
            load(boundaryInfoFilenameSeg, 'boundaryInfo');
            load(sprintf(consts.imageRgbFilename, ii), 'imgRgb');
            load(sprintf(consts.planeDataFilename, ii), 'planeData');
            load(sprintf(consts.watershedFilename, ii), 'pbAll');
            load(sprintf(consts.objectLabelsFilename, ii), 'imgObjectLabels');
            load(sprintf(consts.instanceLabelsFilename, ii), 'imgInstanceLabels');
            load(boundaryFeaturesFilenameSeg, 'boundaryFeatures');

            %[~, instanceLabels] = get_labels_from_instances(boundaryInfo.imgRegions, imgObjectLabels, imgInstanceLabels);

            % Merging
            result = merge_regions(boundaryInfo, boundaryFeatures, classifier, stage, params);
            boundaryInfo = update_boundary_info(boundaryInfo, result, imgRgb, ii, stage);
            save(postMergeFilenameSeg, 'boundaryInfo');

            fprintf(' Done.\n');
        end
    end  
    fprintf('\nFinsished segmenting image #%d.\n', ii);
    fprintf('**************************************************\n');
end

clear global segtestDir;
%clearvars -global segtestDir;
