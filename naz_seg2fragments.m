function segfragments = naz_seg2fragments(boundaryInfo, segmentID)
%NAZ_SEG2FRAGMENTS Obtain list of fragments for given region(s).
%   boundaryInfo - structure output from processBoundaryInfo.m
%   segfragments - cell array, 1-column segmentID, 2-column region points
%% NB! DIDN'T FINISH. SOMETHING IS REALLY WRONG WITH ADJANCY MATRIX IN boundaryInfo.edges.adjacency.

if nargin < 2
    segmentID = unique(boundaryInfo.imgRegions)';
else
    assert(max(segmentID)<=max(boundaryInfo.imgRegions(:)), 'One of IDs exceeing maximum ID number\n');
    assert(length(size(segmentID))==1, '''SegmentID'' must be an array or a single number\n');
end

frag_num = length(boundaryInfo.edges.fragments);
segfragments = cell(length(segmentID), 2);  % 1-column contains region ID, 2-column - region polygon points
c = zeros(1,2);  % polygon centriod
polygon = [];    % polygon points
for i = 1:length(segmentID)
    for k = 1:frag_num
        %DEBUG
        bi = boundaryInfo;
        adj = bi.edges.adjacency{k};
        fprintf(['\n################################\n' ...
                 'k = %d\n%12s: %s(%s)\n%12s: %s\n%12s: %s\n%12s: %s\n'], segmentID,...
                 'adjacency', sprintf('%d ', bi.edges.adjacency{k}), sprintf('%d ',adj),...
                 'superpixs', sprintf('%d ', bi.edges.spLR(k,:)), ...
                 'thetaDir',  sprintf('%.2f', bi.edges.thetaDirected(k)), ...
                 'thetaUndir',sprintf('%.2f', bi.edges.thetaUndirected(k)) ...
               );   
        % if k-th fragment is adjacent to i-th segment (from any side - left or right)
        if any(boundaryInfo.edges.spLR(k,:)==segmentID(i));
            % since don't know in advance number of fragments per segment,
            % each cell is populated dynamically, without size pre-allocation
            polygon = [polygon; boundaryInfo.edges.fragments{1,k}];

            %DEBUG
            fprintf('fragments:\n');
            frag = boundaryInfo.edges.fragments{1,k};
            disp(frag);
            
            if ~isempty(adj)
                [c r] = ind2sub(size(bi.imgRegions),adj');
                ha1 = scatter(r(1), c(1), 'og');
                ha2 = scatter(r(2:end), c(2:end), 'or');
            end
            
            hf = plot(frag(:,1), frag(:,2), 'y', 'LineWidth', 0.5);
            delete([hf ha1 ha2]);
        end
%     clc;
    end

segfragments{i,1} = segmentID(i);
segfragments{i,2} = polygon;

hd = plot(polygon(:,1), polygon(:,2), 'y', 'LineWidth', 1);
delete(hd);

polygon = [];
end

end

% boundaryInfo structure does not include fragments on the image border.
% These borders are needed to close the polygon (for inpolygon function)
function bfrags = border_fragments(imgRegions)
    % imgRegions is expected to be passed from boundaryInfo.imgRegions==k,
    % where k is needed segmentID
    

end

% NB! Didn't work!
% sorting polygon points clockwise to obtain correct poins for plot() function.
% % the centroid:
% c(1) = mean(polygon(:,1));
% c(2) = mean(polygon(:,2));
% % angle of each point:
% 
% th = atan2(polygon(:,2) - c(2), polygon(:,1) - c(1));
% % correct sort order:
% [~, ndxs] = sort(th);
% poly_x = polygon(:,1);
% poly_y = polygon(:,2);
% polygon(:,1) = poly_x(ndxs);
% polygon(:,2) = poly_y(ndxs);
