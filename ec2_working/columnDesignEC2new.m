function [h, b, noRebar, phiRebar, areaRebar, M_Rd, phiShear, spacingShear] = columnDesignEC2new(fck, fyk , cover, N_Axial, My_h, Mz_b, width)

%data = 'dataColumns.csv' ; %input('Filename: ');
%unsData = importdata(data);
abaco = importdata('info\abacusC12_50S500quarter.mat');
longReinforce = importdata('info\steel_column.csv');
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5; 
%fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
%fywd = fyd;
%% EC 2 Design
%longitudinal reinforcement

if nargin == 7
    h = width;
    b = width;
else
    h = .2;
    b = .2;
end

diffe = -1;
while diffe < 0
    redAxial = N_Axial / (b * h * fcd * 1000);
    redBenMom1 = Mz_b / (b * h^2 * fcd * 1000);
    redBenMom2 = My_h / (h * b^2 * fcd * 1000);
    
    reinfPerc = abaco(redAxial, redBenMom1, redBenMom2);
    AsAbacus = reinfPerc * b * h * fcd / fyd;
    
    Asmin = max([.1 * N_Axial / (fyd * 1000), .002 * b * h]);
    reinfArea = max([Asmin, AsAbacus]);
  
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
    areaRebar = longReinforce(minIndex,3);
    
    diffe = b - 2 * (cover + .02) - (phiRebar/1000 * (noRebar/4 + 1)) - (max([.02, phiRebar/1000]) * (noRebar/4));
    b = b + .05;
    h = h + .05;
end
b = b - .05;
h = h - .05;

%% M_Rd 
reinfPercFin = areaRebar * fyd / (b * h * fcd);
redAxialFin = N_Axial / (b * h * fcd * 1000);
diff = [];
for i = 0 : .005 : 1
    reinfPerc = abaco(redAxialFin, i, i);
    diff = [diff; [i , abs(reinfPercFin - reinfPerc)]];
end

[value1, index1] = min(diff(:, 1));
redBenMom = (index1 - 1) * .005;
M_Rd = redBenMom * h * b^2 * fcd * 1000 ; 

%% stirrups
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