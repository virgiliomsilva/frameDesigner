function [h, b, noRebar, phiRebar, areaRebar, M_Rd] = columnDesignDC3(fck, fyk , cover, N_Axial, My_h, Mz_b, b_input, h_input, phi)
abaco = importdata('info\abacusC12_50S500quarter.mat');
longReinforce = importdata('info\steel_columnEC8.csv');

fcd = fck / 1.5; 
fyd = fyk / 1.15;

%% longitudinal reinforcement

if nargin >= 8
    b = ceil(b_input * 20) / 20;
    h = ceil(h_input * 20) / 20;
else
    b = .2;
    h = .2;
end
    
redAxial = N_Axial / (b * h * fcd * 1000);

while redAxial > .55
    b = b + .05;
    h = h + .05;
    redAxial = N_Axial / (b * h * fcd * 1000);
end

diff = -1;
while diff < 0
    redAxial = N_Axial / (b * h * fcd * 1000);
    redBenMom1 = Mz_b / (b * h^2 * fcd * 1000);
    redBenMom2 = My_h / (h * b^2 * fcd * 1000);
    
    reinfPerc = abaco(redAxial, redBenMom1, redBenMom2);
    AsAbacus = reinfPerc * b * h * fcd / fyd;
    
    Asmin = max([.1 * N_Axial / (fyd * 1000), .002 * b * h]);
    reinfArea = max([Asmin, AsAbacus]);
    
    if nargin == 9
        longReinforce(longReinforce(:,1) <= phi, :) = [];
    end

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
    
    diff = b - 2 * (cover + .02) - (phiRebar/1000 * (noRebar/4 + 1)) - (max([.02, phiRebar/1000]) * (noRebar/4));
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
    diff = [diff; [i, abs(reinfPercFin - reinfPerc)]];
end

[value index] = min(diff(:,2));
redBenMom = (index - 1) * .005;
M_Rd = redBenMom * h * b^2 * fcd * 1000 ; 