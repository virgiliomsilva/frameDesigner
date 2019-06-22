%% MINIMUM WIDTH FINDER
% Given a bar (column) name will find the minimum width possible compatible 
% with the adjacent bars (beams)
function [minWidth] = minWidFind(barName, element, beams)
nodes = element(element(:,1) == barName,[2 3]);
bars = [];
for i = 1 : 2
    for j = 1 : size(element,1)
        if (element(j, i+1) == nodes(1) | element(j, i+1) == nodes(2)) & element(j,4) ~= 3
            bars(end+1) = element(j,1);
        end
    end 
end
bars = unique(bars);
widths = [];
for i = 1 : length(bars)
    widths(end+1, 1) = beams(beams(:,1) == bars(i) , 3);
end
minWidth = min(widths);
end