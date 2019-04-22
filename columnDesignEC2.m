function [b, noRebar, phiRebar, phiShear, spacingShear] = columnDesignEC2(fck, fyk , cover, N_Axial, My_h, Mz_b, b_input)

%data = 'dataColumns.csv' ; %input('Filename: ');
%unsData = importdata(data);
abaco = importdata('info\abacusC12_50S500A1_composta.csv');
longReinforce = importdata('info\steel_column.csv');
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5; 
%fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
%fywd = fyd;
%% EC 2 Design
%longitudinal reinforcement

if nargin == 6
    b = max(.2, round(sqrt(N_Axial / (fcd * 1000)) * 20)/20);
else
    b = b_input
end

diff = -1;

while diff < 0
    redAxial = N_Axial / ( b^2 * fcd * 1000);
    benMom = sqrt(My_h^2 + Mz_b^2);
    redBenMom = benMom / (b^3 * fcd * 1000);
    for j = 1 : size(abaco,2)
        mAux(1,j) = abaco(1,j) - redAxial;
    end
    [value col]= min(abs(mAux));
    for k = 1 : size(abaco,1)
        mAux2(k,1) = abaco(k,1) - redBenMom;
    end
    [value row] = min(abs(mAux2));
    reinfPerc = abaco(row, col);
    Asmin = max([.1 * N_Axial / (fyd * 1000), .002 * b^2]);
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

sAvaila = unique(shearReinforce(:,3));
spacing = min([15 * phiRebar/1000, b, .3]);
phiShear = 8;

for k = 1 : size(sAvaila,1)
    if sAvaila(k) - spacing <= 0
        auxMatrix(k,1) = sAvaila(k) - spacing;
    else
        auxMatrix(k,1) = -1000;
    end
end

[value index] = max(auxMatrix);
spacingShear = sAvaila(index);

%final data
%reinfSolu(i,1) = finalData(i,1);
% reinfSolu(i,2) = b;
% reinfSolu(i,3) = noRebar;
% reinfSolu(i,4) = phiRebar;
% 
% reinfSolu(i,5) = 8;
% reinfSolu(i,6) = spacingShear;

%reinfSolu(i,8) = longReinforce(minIndex,3) * (b - 2 * (cover + .02)) * fyd * 1000

%reinfSolu(i,8) = longReinforce(minIndex,3) * fyd / (b * b * fcd);
%  reinfSolu(i,9) = N_Axial / (b * b * fcd * 1000);
% reinfSolu(i, 10) = My_h / ( b * b * b * fcd * 1000);
% reinfSolu(i, 11) = Mz_b / ( b * b * b * fcd * 1000);


       