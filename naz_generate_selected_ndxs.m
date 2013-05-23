function indices = naz_generate_selected_ndxs(match_str)
%NAZ_GENERATE_NDXS Generate indices based on match string (with wildcards)
% match_str expects to match 6-digit string, e.g.: "2*", "*001*", "10000?".

Consts;
glob = strcat(consts.imageRgbDir, 'rgb_', match_str, '.png');
filelist = dir(glob);

indices = {filelist.name}';
i = regexp(indices{1}, '\d{6}');
indices = cellfun(@(x) str2double(x(i:i+5)), indices);
indices = sprintf('[%s ]', sprintf(' %06d', indices));

end

