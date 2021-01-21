function [id] = get_id(i,j,str)
% i and j are row and collumn of main node. str is optional. without it you
% get the main node id. set it to 'L','R','U'' or 'D' to get a secondary
% node id

if nargin<3
    id = strcat(num2str(i), num2str(j));
else
    id = strcat(num2str(i), num2str(j), str);
end


end

