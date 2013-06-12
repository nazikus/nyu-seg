aa.pair_ind = bi.edges.spLR;
aa.coord = bi.junctions.position;
aa.regions = bi.imgRegions;
aa.centroids = regionprops(aa.regions, 'Centroid');

set(0,'DefaultFigureWindowStyle','docked');
figure, imshow(aa.regions,[])%, colormap jet;
hold on
for k = 1:numel(aa.centroids)
    c = aa.centroids(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'Bold', 'Color', [0 0.8 0]);
end

% for k=1:numel(bi.edges.fragments);
%     frag = bi.edges.fragments{k};
%     plot(frag(:,1), frag(:,2), 'y', 'LineWidth', 0.5);
%     hline = plot(frag(:,1), frag(:,2), 'r', 'LineWidth', 3);
%     fprintf(['\n################################\n' ...
%                'k = %d\n' ...
%                '%12s: %s\n', ...
%                '%12s: %s\n', ...
%                '%12s: %s\n', ...
%                '%12s: %s\n', ...
%             ], ...
%                'adjacency', sprintf('%d ', bi.edges.adjacency{k}), ...
%                'superpixs', sprintf('%d ', bi.edges.spLR(k,:)), ...
%                'thetaDir',  sprintf('%.2f', bi.edges.thetaDirected(k)), ...
%                'thetaUndir',sprintf('%.2f', bi.edges.thetaUndirected(k)) ...
%            );   
%                
%     delete(hline);
%     clc;
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size(aa.regions);
[ry rx] = ind2sub(sz, randsample(prod(sz),150));

for r = unique(aa.regions)';
    fprintf('region = %d\n', r);
    %---------------
    % padding & removing outliers, to avoid odd behavior of 'contour' function (bug?)
    reg = padarray(aa.regions==r, [1 1]);
    cont_ini = contourc(double(reg),1);
    [~, cont_outlier_c] = find(cont_ini<1);
    %[~, cont_outlier_c] = ind2sub(size(cont_ini), cont_outlier);
    cont_c = logical(1:length(cont_ini));
    cont_c(cont_outlier_c) = false;
    cont = cont_ini(:,cont_c) - 1; % -1 to compensate coordinates shift because of padding
    clear cont_outlier cont_outlier_c cont_c;
    %---------------
    
    in = inpolygon(rx,ry,cont(1,:),cont(2,:));
    
    hc = plot(cont(1,:), cont(2,:), 'r');
    hr = plot(rx(in), ry(in), 'mo', rx(~in), ry(~in), 'cx');
    if any(r==[32 62])
        
    end
    delete(hr);
    delete(hc);
    
end
hold off

% NB! Didn't work, polygon is not convex (squared edges)
% -------
% for k = unique(aa.regions)'
%    bw = im2bw(aa.regions == k);
%    bnd = bwboundaries(bw);
%    plot(poly(:,1), poly(:,2), 'y', 'LineWidth', 1);
%    pause;
% end

% NB! DIDN'T FINISH. SOMETHING IS REALLY WRONG WITH ADJANCY MATRIX IN boundaryInfo.edges.adjacency.
% segfrags = naz_seg2fragments(bi);
% for m = 1:length(segfrags)
%     frag = segfrags{m,2};
%     plot(frag(:,1), frag(:,2), 'y', 'LineWidth', 0.5);
%     hseg = plot(frag(:,1), frag(:,2), 'r', 'LineWidth', 2);
%     fprintf('\n################################\nk = %d\n', segfrags{m,1} );   
%     
%     delete(hseg);
%     clc;
% end
