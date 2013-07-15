function boundaryOutput = naz_remove_fragments(boundaryInput, fragmentList)
%NAZ_REMOVE_FRAGMENTS Remove fragments and merges corresponding regions
% NB! Reliably alters only boundaryInfo.edges.fragments (used in naz_plot_fragments),
% and boundaryInfo.imgRegions (used in evaluate_segmentation); If you need to use other information 
% (like junctions) from penalized boundaryInfo, you need to remove it properly in this funciton.

boundaryOutput = boundaryInput;
% removing image region (labels) by assigning region 
% on the left from fragment label of the region on the right.
for k=1:size(fragmentList,1)
    f = fragmentList{k,1};
    for p = 1:size(f,1)
        L = f{p}(1);
        ndxs = boundaryOutput.imgRegions == L;
        boundaryOutput.imgRegions(ndxs) = f{p}(2);
        % changing paired labels ahead of fragmentlist accordingly
        for v=p:size(f,1)
            if any(f{v}==L)
                f{v}( f{v}==L ) = f{p}(2);
            end            
        end
    end
end

% after the previous loop, we have labeling "gaps" (e.g. 1 2 3 . 6 7 8)
% here we will shift (redefine) all the labeling in sequential order
labels = unique(boundaryOutput.imgRegions);
for k = 1:length(labels)
    boundaryOutput.imgRegions(boundaryOutput.imgRegions==labels(k)) = k;    
end
clear lables;

% removing fragments
frags = []; %flatenning cells in single-column array
for k = 1:size(fragmentList,1)
    f = fragmentList{k,2};
    for p = 1:size(f,1)
        frags = [frags; f{p}(:)]; %#ok<AGROW>
    end
end
frags = sort(frags); 
ndxs = true(size(boundaryInput.edges.fragments));
ndxs(frags) = false;
boundaryOutput.edges.fragments = boundaryInput.edges.fragments(ndxs);

% removing some other fields (with no purpose so far, just for consistency)
%boundaryOutput.ne = boundaryOutput.ne - length(frags);
%boundaryOutput.edges.indices   = boundaryInput.edges.indices(ndxs);
%boundaryOutput.edges.spLR      = boundaryInput.edges.spLR(ndxs',:);
%boundaryOutput.edges.thetaDirected   = boundaryInput.edges.thetaDirected(ndxs');
%boundaryOutput.edges.thetaUndirected = boundaryInput.edges.thetaUndirected(ndxs');

% NB! need to take care of corresponding boundaryInfo.junctions, before removing boundaryInfo.edges.juncitons
%boundaryOutput.edges.junctions = boundaryInput.edges.junctions(ndxs');
% NB! length(adjacency) == ne*2
%boundaryOutput.edges.adjacency = boundaryInput.edges.adjacency(ndxs'); 
end

