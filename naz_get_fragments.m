function [spPairNdxs fragmentNdxs] = naz_get_fragments(regionarr, boundaryInfo)
%NAZ_GET_FRAGMENTS Obtain list of fragments indices between specified regions (superpixels)

    sp = sort(boundaryInfo.edges.spLR, 2); % sort for comparison (we don't care about left-right relation in here)
    sp_pairs = sort( combnk(regionarr,2), 2);
    
    fragmentNdxs = cell(size(sp_pairs,1), 1);
    spPairNdxs     = cell(size(sp_pairs,1), 1);
    for i = 1:size(sp_pairs, 1);
         ndxs = repmat(sp_pairs(i,:), [size(sp,1) 1]) == sp;  % find logical indices of fragments that are sepearting sp-regions
         ndxs = find(ndxs(:,1) & ndxs(:,2)); % convert logical indices to numerical
         
         spPairNdxs{i}     = sp_pairs(i,:);
         fragmentNdxs{i} = ndxs;
    end
end

