function newMatch = naz_dissimilar_region_frags(fragmentList, LR, rm, threshold)
fragList = fragmentList{LR};
newMatch = fragList;
for k = 1:size(fragList,1)
    f = fragList{k,1};
    c = 0;
    for p = 1:size(f,1)
        % region labels
        rR = fragList{k,3};
        rLl = f{p}(1); %region in the left image on the left  side of the current fragment
        rLr = f{p}(2); %region in the left image on the right side of the current fragment
        
        % sift indices
        siR = rm.region2ind{3-LR}(rR,:);  % '3-LR' means if LR==1, then its 2, if LR==2, then its 1.
        siLl = rm.region2ind{LR}(rLl,:);
        siLr = rm.region2ind{LR}(rLr,:);
        
        % sift counts
        scLlR = sum(siLl&siR);
        scLl = sum(siLl);
        scLrR = sum(siLr&siR);
        scLr = sum(siLr);
        
        % region similarity between one of the left region adn the right region
        simLl = scLlR/scLl;
        simLr = scLrR/scLr;
        if simLl > threshold && simLr > threshold
            fprintf('Similarity found! rR: %d (#%d), frag: %d (Ll:%d=Lr:%d)\n', rR, k, p, rLl, rLr);
        else        
            newMatch{k,1}(p-c,:) = [];
            newMatch{k,2}(p-c,:) = [];
            c = c + 1;
        end
    end    
end

end