function cont = naz_region_contour(bwregion)
%NAZ_REGION_CONTOUR Fixed matlab contour function
    %FIX some awkward bug, incorrect contour sometimes
    
    % padding & removing outliers, to supress odd results from 'contour' function (bug?)
%     reg = padarray(bwregion, [1 1]);
    [cont labels] = bwboundaries(bwregion, 8, 'noholes');
    lens = cellfun(@length, cont);
    [~, ind] = max(lens);
    cont = flipud(cont{ind}');
end

% hh = plot(cont_ini(1,:)+429, cont_ini(2,:), 'Color', 'c', 'LineWidth', 1);
% aa = bwboundaries(bwregion, 8);




%     % padding & removing outliers, to supress odd results from 'contour' function (bug?)
%     reg = padarray(bwregion, [1 1]);
%     cont_ini = contourc(double(reg),1);   % inital values obtain from 'contour' func, with odd points <1
%     [~, cont_outlier_c] = find(cont_ini<1);
%     cont_c = logical(1:length(cont_ini));
%     cont_c(cont_outlier_c) = false;
%     cont = cont_ini(:,cont_c) - 1; % -1 to compensate coordinates shift because of padding
