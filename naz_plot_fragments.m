function naz_plot_fragments(fragNdxs, boundaryInfo, shift, color, width)
%NAZ_PLOT_EDGES Plot given edges in the current figure with paired images, hence the 'shift'

    conf.sizeY = 0;
    conf.imgGap = shift; % <-- lame workaround
    
    frags = [];
    for j = 1:length(fragNdxs)
        frags = [frags; fragNdxs{j}(:)]; %#ok<AGROW>
    end
    f = boundaryInfo.edges.fragments( frags );  % get fragmentlist
    naz_plot_paired_edges({f, []}, conf, color, width);

end

