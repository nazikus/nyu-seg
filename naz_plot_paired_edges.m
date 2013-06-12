function naz_plot_paired_edges(edgesIm, conf, color, width, orientation)
%NAZ_PLOT_EDGES Plot segmentation edges in current figure
switch (nargin)
    case 2
        color = 'r';
        width = 0.5;
        % orientation = ?
    case 3
        width = 0.5;
        % orientation = ?
    case 4
        % orientation = ?
end

imgNum = length(edgesIm);
for i = 1:imgNum
   if i==1
       shift = 0;
   else
       shift = conf.sizeY + conf.imgGap;
   end
   edges = edgesIm{i};
   for k = 1:length(edges)
        plot(edges{k}(:,1)+shift, edges{k}(:,2), 'Color', color, 'LineWidth', width);
   end
end
end

