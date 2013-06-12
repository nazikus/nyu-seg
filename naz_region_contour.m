function cont = naz_region_contour(bwregion)
%NAZ_REGION_CONTOUR Fixed matlab contour function
    %FIX some awkward bug, incorrect contour sometimes
    
    % padding & removing outliers, to supress odd results from 'contour' function (bug?)
    reg = padarray(bwregion, [1 1]);
    cont_ini = contourc(double(reg),1);   % inital values botain from 'contour' func, with odd points <1
    [~, cont_outlier_c] = find(cont_ini<1);
    cont_c = logical(1:length(cont_ini));
    cont_c(cont_outlier_c) = false;
    cont = cont_ini(:,cont_c) - 1; % -1 to compensate coordinates shift because of padding
end

