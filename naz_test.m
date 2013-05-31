aa.pair_ind = bi.edges.spLR;
aa.coord = bi.junctions.position;
aa.regions = bi.imgRegions;

aa.centroids = regionprops(aa.regions, 'Centroid');
figure, imshow(aa.regions,[]);
hold on
for k = 1:numel(aa.centroids)
    c = aa.centroids(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'Bold', 'Color', [0.7 0 0]);
end

for k=1:numel(bi.edges.fragments);
    frag = bi.edges.fragments{k};
    plot(frag(:,1), frag(:,2), 'y', 'LineWdith', 0.5);
    fprintf(['\n################################\n' ...
               'k = %d\n' ...
               '%12s: %s\n', ...
               '%12s: %s\n', ...
               '%12s: %s\n', ...
               '%12s: %s\n', ...
            ], ...
               'adjacency', sprintf('%d ', bi.edges.adjacency{k}), ...
               'superpixs', sprintf('%d ', bi.edges.spLR(k,:)), ...
               'thetaDir',  sprintf('%.2f', bi.edges.thetaDirected(k)), ...
               'thetaUndir',sprintf('%.2f', bi.edges.thetaUndirected(k)) ...
           );   
               
    pause;
    clc;
end

hold off
clear k c poly bw bnd frag;



% for k = unique(aa.regions)'
%    bw = im2bw(aa.regions == k);
%    bnd = bwboundaries(bw);
%    plot(poly(:,1), poly(:,2), 'y', 'LineWidth', 1);
%    pause;
% end


