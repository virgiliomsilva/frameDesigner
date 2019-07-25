function [longReinforce] = columnRebarMin(areaMin)
longReinforce = importdata('info\steel_columnEC8.csv');
longReinforce = longReinforce(:,[1:3]);
longReinforce = longReinforce(longReinforce(:,3) >= areaMin, :)
end