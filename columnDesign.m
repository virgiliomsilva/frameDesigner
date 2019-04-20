 %COLUMN DESIGN
clear
clc

%% import data and inputs


data = 'dataColumns.csv' ; %input('Filename: ');
unsData = importdata(data);
tic
abaco = importdata('abacusC12_50S500A1_composta.csv');
longReinforce = importdata('steel_column.csv');
%shearReinforce = importdata('steel_shear.csv');

fck = 25;  %MPa
fyk = 500;  %MPa
cover = .035; 

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd;
%% transform data
[noRow, noCol] = size(unsData.textdata);
barNameAux = extractBefore(unsData.textdata(2:noRow,:),"/");

for i = 1 : (noRow - 1)
    auxMat(i, 1) = str2double(barNameAux(i, 1));
end

sortData = abs([auxMat, unsData.data]);

row = 1;
rowz = 1;
for j = unique(sortData(:, 1)).'
    for i = 1 : (noRow - 1) 
        if j == sortData(i, 1)
        auxBlock(row, :) = sortData(i, :);
        row = row + 1;
        end
    end
row = 1;
auxMax = max(auxBlock, [], 1);
finalData(rowz, :) = auxMax;
rowz = rowz + 1;
end

noBars = rowz - 1 ;
%reinfSolu = zeros(noBars, ###) %%%%%%%%%%%%%%%%%%%%%
clear rowNum colNum unsData barNameAux auxMat sortData row auxBlock auxMax...
    i j rowz data noCol noRow
%% EC 2 Design

for i = 1 : noBars
    
    %longitudinal reinforcement
    b = max(.2, round(sqrt(finalData(i,2) / (fcd * 1000)) * 20)/20);
    diff = -1;
    while diff < 0
        redAxial = finalData(i, 2) / ( b^2 * fcd * 1000);
        benMom = sqrt(finalData(i, 5)^2 + finalData(i, 6)^2);
        redBenMom = benMom / (b^3 * fcd * 1000);
        for j = 1 : 13
            mAux(1,j) = abaco(1,j) - redAxial;
        end
        [value col]= min(abs(mAux));
        for k = 1 : 101
            mAux2(k,1) = abaco(k,1) - redBenMom;
        end
        [value row] = min(abs(mAux2));
        reinfPerc = abaco(row, col);
        Asmin = max([.1 * finalData(i,2) / (fyd * 1000), .002 * b^2]);
        AsAbacus = reinfPerc * b^2 * fcd / fyd;
        reinfArea = max(Asmin, AsAbacus * 2);
    
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfArea > 0
                diffAux(j,1) = longReinforce(j,3) - reinfArea;
            else
                diffAux(j,1) = 1000;
            end
        end
    
        [minDiff, minIndex] = min(diffAux);
        noRebar = longReinforce(minIndex,2);
        phiRebar = longReinforce(minIndex,1);
    
        diff = b - 2 * (cover + .02) - (phiRebar/1000 * (noRebar/4)) - (max([.02, phiRebar/1000]) * (noRebar/4 - 1));
        b = b + .05;
    end
    b = b - .05;
    
    %stirrups
    %phi 8mm 
    
    spacing = min([15 * phiRebar/1000, b, .3]);
    
    for k = 1 : 7
        if sAvaila(k) - spacing <= 0
            auxMatrix(k,1) = sAvaila(k) - spacing;
        else
            auxMatrix(k,1) = -1000;
        end
    end
    
    [value index] = max(auxMatrix);
    
    %final data
    reinfSolu(i,1) = finalData(i,1);
    reinfSolu(i,2) = b;
    reinfSolu(i,3) = noRebar;
    reinfSolu(i,4) = phiRebar;
    
    reinfSolu(i,5) = 8;
    reinfSolu(i,6) = sAvaila(index);
    
    reinfSolu(i,8) = longReinforce(minIndex,3) * (b - 2 * (cover + .02)) * fyd * 1000
    
    %reinfSolu(i,8) = longReinforce(minIndex,3) * fyd / (b * b * fcd);
  %  reinfSolu(i,9) = finalData(i,2) / (b * b * fcd * 1000);
   % reinfSolu(i, 10) = finalData(i,5) / ( b * b * b * fcd * 1000);
   % reinfSolu(i, 11) = finalData(i,6) / ( b * b * b * fcd * 1000);
end


%% EC 8 design

%   sAvaila = unique(shearReinforce(:,3));

for i = 1 : noBars
    
    %DC3
    %longitudinal reinforcement
    b = max([.25, round(sqrt(finalData(i,2) / (fcd * 1000)) * 20)/20, round(sqrt(finalData(i,2)/(.55 * fcd * 1000))) * 20/20]); 
    diff = -1;
    while diff < 0
        roL = 1;
        while roL > .04
            redAxial = finalData(i, 2) / ( b^2 * fcd * 1000);
            benMom = sqrt(finalData(i, 5)^2 + finalData(i, 6)^2);
            redBenMom = benMom / (b^3 * fcd * 1000);
            for j = 1 : 13
                mAux(1,j) = abaco(1,j) - redAxial;
            end
            [value col]= min(abs(mAux));
            for k = 1 : 101
                mAux2(k,1) = abaco(k,1) - redBenMom;
            end
            [value row] = min(abs(mAux2));
            reinfPerc = abaco(row, col);
            Asmin = max([.1 * finalData(i,2) / (fyd * 1000), .002 * b^2]);
            AsAbacus = reinfPerc * b^2 * fcd / fyd;
            reinfArea = max(Asmin, AsAbacus * 2);

            for j = 1 : size(longReinforce,1)
                if longReinforce(j,3) - reinfArea > 0 & longReinforce(j,2) >=8
                    diffAux(j,1) = longReinforce(j,3) - reinfArea;
                else
                    diffAux(j,1) = 1000;
                end
            end

            [minDiff, minIndex] = min(diffAux);
            noRebar = longReinforce(minIndex,2);
            phiRebar = longReinforce(minIndex,1);
            areaRebar = longReinforce(minIndex,3);

            diff = b - 2 * (cover + .02) - (phiRebar/1000 * (noRebar/4)) - (max([.02, phiRebar/1000]) * (noRebar/4 - 1));
            b = b + .05;
            roL = areaRebar / b^2;
        end
    end
    b = b - .05;
    
    reinfSolu(i,10) = b;
    reinfSolu(i,11) = noRebar;
    reinfSolu(i,12) = phiRebar;
    
    
    
    %{
    
    
    %stirrups
    %phi 8mm 
    
    spacing = min([15 * phiRebar/1000, b, .3]);
    
    for k = 1 : 7
        if sAvaila(k) - spacing <= 0
            auxMatrix(k,1) = sAvaila(k) - spacing;
        else
            auxMatrix(k,1) = -1000;
        end
    end
    
    [value index] = max(auxMatrix);
    
    %final data
    reinfSolu(i,1) = finalData(i,1);
    reinfSolu(i,2) = b;
    reinfSolu(i,3) = noRebar;
    reinfSolu(i,4) = phiRebar;
    
    reinfSolu(i,5) = 8;
    reinfSolu(i,6) = sAvaila(index);
    
    reinfSolu(i,8) = longReinforce(minIndex,3) * (b - 2 * (cover + .02)) * fyd * 1000
    
    %reinfSolu(i,8) = longReinforce(minIndex,3) * fyd / (b * b * fcd);
  %  reinfSolu(i,9) = finalData(i,2) / (b * b * fcd * 1000);
   % reinfSolu(i, 10) = finalData(i,5) / ( b * b * b * fcd * 1000);
   % reinfSolu(i, 11) = finalData(i,6) / ( b * b * b * fcd * 1000);
    %}
end



toc
            