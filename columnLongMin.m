%%
% if a second argument is inserted, will also give patterns with only 4 rebars
function [longReinforce] = columnLongMin(areaMin, flag)
    if nargin == 1
        longReinforce = importdata('info\steel_columnEC8.csv');
        longReinforce = longReinforce(:,[1:3]);
        longReinforce = longReinforce(longReinforce(:,3) >= areaMin, :);
    elseif nargin == 2
        longReinforce = importdata('info\steel_columnDC1.csv');
        longReinforce = longReinforce(:,[1:3]);
        longReinforce = longReinforce(longReinforce(:,3) >= areaMin, :);
    end
end