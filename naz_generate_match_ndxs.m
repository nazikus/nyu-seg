function matchNdxs = naz_generate_match_ndxs( group_arr )
%NAZ_GENERATE_MATCHNDXS Generate indices pairs from 'consts.useNdx' according to grouping.
% group_arr contains indices that correspond to range end. 
% E.g.: if we have a range 001 - 015, 016 - 020, 020 - 100, then
% 'group_arr' will contain = [015 020 100];

%% SIMPLIFIED VERSION, CREATE ALL PAIRS IN PAIRWISE MANNER (ignoring 'group_arr' argument)
Consts;
Ndxs = consts.useNdx;
matchNdxs = cell(length(Ndxs)-1, 1);
for i = 2:length(Ndxs)
    matchNdxs{i-1} = [Ndxs(i-1) Ndxs(i)];
end

%%%%%%%%%%%%%%% EXLCLUDE PAIRS CROSSING THE GROUP RANGES %%%%%%
% Consts;
% Ndxs = consts.useNdx;
% matchNdxs = cell(length(Ndxs) - length(group_arr), 1);
% k = 1;
% for i = 2:length(Ndxs)
%     if Ndxs(i) <= group_arr(1)
%         matchNdxs{k} = [Ndxs(i-1) Ndxs(i)];
%         k = k + 1;
%     else
%        group_arr = group_arr(2:end); 
%     end
% end

end

