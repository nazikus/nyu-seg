% INDOOR_SCENE_SEG_SUP
%
% Files
%   Consts                                          - Contains the constants used throughout the project.
%   Params                                          - Contains the parameters to the different modules of the
%   README                                          - Requirements:
%   run_compile_all                                 - 
%   run_create_dataset_rgbd_sift                    - Creates a single (giant) dataset of all the SIFT descriptors extracted from the dataset. The
%   run_create_dataset_structure_class_features_gt  - Creates a dataset comprised of structure class features extracted from the ground truth (manually
%   run_create_dataset_structure_class_features_seg - Creates a dataset comprised of structure class features extracted from the segmented (bottom-up) 
%   run_create_dataset_support_features_gt          - Creates a dataset of support features to use for training a support
%   run_evaluate_segmentation                       - Evaluates the segmentation pipeline using the overlap metric from Hoiem et al's Recovering Occlusion Boundaries
%   run_evaluate_support                            - 
%   run_extract_regions_from_labels                 - Create the regions from the ground truth (manually annotated) labels.
%   run_extract_sift_features                       - Extract SIFT features from each frame's RGB and Depth image.
%   run_extract_structure_class_features_gt         - Extracts and saves all of the region-to-structure class features using the ground truth (manually
%   run_extract_structure_class_features_seg        - Extracts and saves all of the region-to-structure class features.
%   run_extract_support_features_gt                 - Extracts support features for every image in the dataset using ground
%   run_extract_support_features_seg                - Extracts support features for every image in the dataset using ground
%   run_extract_surface_normals                     - Loads all of the in-painted depth images and estimates the surface
%   run_label_watershed_segments                    - 
%   run_learn_sc_dict_from_sift                     - Learns a dictionary of atoms used for sparse coding SIFT descriptors.
%   run_pipeline_segmentation                       - Contains the scripts that implement the segmentation pipeline.
%   run_pipeline_support_inference_gt               - Runs the segmentation pipeline followed by support inference on the
%   run_pipeline_support_inference_seg              - Runs the segmentation pipeline followed by support inference on the
%   run_regions2labels                              - Determines the class and instance labels for each of the initial superpixel segments produced by
%   run_rgbd2planes                                 - Finds major scene surfaces and rotates the entire scene such that the
%   run_save_individual_files                       - Loads the labeled dataset and saves it out to individual files which will
%   run_show_images                                 - Visualizes all of the images.
%   run_support_inference_image_plane_rules_gt      - Evaluates the accuracy of Support Baseline #1, a rule based approach
%   run_support_inference_image_plane_rules_seg     - Evaluates the accuracy of Support Baseline #1, a rule based approach
%   run_support_inference_ip_gt                     - Performs inference by finding the optiminal binary (IP) solutions on ground truth (human
%   run_support_inference_lp_gt                     - Infers the supporting regions, support directions and structure classes
%   run_support_inference_lp_seg                    - Infers the supporting regions, support directions and structure classes
%   run_support_inference_structure_class_rules_gt  - Evaluates the accuracy of Support Baseline #2, a set of rules based on
%   run_support_inference_structure_class_rules_seg - Evaluates the accuracy of Support Baseline #2, a set of rules based on
%   run_support_inference_support_classifier_gt     - Evaluates the accuracy of Support Baseline #3, an approach in which the
%   run_support_inference_support_classifier_seg    - Evaluates the accuracy of Support Baseline #3, an approach in which the
%   run_train_boundary_classifiers                  - Trains several stages of boundary classifiers
%   run_train_boundary_classifiers_with_support     - Trains several stages of boundary classifiers.
%   run_train_floor_classifier_gt                   - Trains a classifier to predict the dominant label of the classifier.
%   run_train_floor_classifier_seg                  - Trains a classifier to predict the dominant label of the classifier.
%   run_train_structure_class_classifier_gt         - Trains a classifier to predict a region's Structure Class using the features extracted from the
%   run_train_structure_class_classifier_seg        - Trains a classifier to predict a region's Structure Class using the features extracted from the
%   run_train_support_classifier_gt                 - 
%   run_watershed_segmentation                      - Performs an initial watershed segmentation on each RGBD image.
