%% COMPATIBILITY FUNCTION
% given a steel reinforcement pattern (de facto area) the function returns
% the other possible patterns suitable to merge with the given on
% flag inputs may be 'EC2' or 'EC8'
function [newComp] = columnCompEC8(token, flag)
% token = 0.003216991 ;
% flag = 'Ec8' ;
table = importdata('info\columnCompatibility.csv');
%flag c
% if strcmpi(flag, 'EC8')
%     table = importdata('info\steel_columnEC8.csv');
% elseif strcmpi(flag, 'EC2')
%     table = importdata('info\steel_column.csv');
% end
%possible diameters
% phiTable = [12, 16, 20, 25, 32];

% [~, col] = find(table(:,3) == token);
[~, col] = find(table == token);
newComp = table(:, [col-2 : col]);


