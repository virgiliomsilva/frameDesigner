%% hardest way to get the MRD of a beam!!
function [MRD] = MrdBeam (fck, fyk, cover, b, h, A1, A2) 
if (fyk == 500 && fck > 12 && fck < 50)
    abaco = importdata('info\abacus_C12_50A500_allBend_int.mat');
elseif (fyk == 400 && fck > 12 && fck < 50)
    abaco = importdata('info\abacus_C12_50A400_allBend_int.mat');
else
    error('Materials pair not supported!')
end

fyd = fyk / 1.15;
fcd = fck / 1.5;

Aratio = min([A1, A2]) / max([A1, A2]);

d = h - (cover + .02);
reinfPerc = max([A1, A2]) * fyd / (b * d * fcd);

Ind1 = find(abaco(1,:) < Aratio, 1, 'last'); 
ratio1 = abaco(1, Ind1);
Ind2 = find(abaco(1,:) >= Aratio, 1); 
ratio2 = abaco(1, Ind2);

for i = 2 : size(abaco,1)
     auxR1(i-1) = abs(abaco(i, Ind1) - reinfPerc);
     auxR2(i-1) = abs(abaco(i, Ind2) - reinfPerc);
end

[~, min1] = min(auxR1);
redBenMom1 = abaco(min1 + 1, 1);
[~, min2] = min(auxR2);
redBenMom2 = abaco(min2 + 1, 1);

redBenMom = interp1([ratio1, ratio2], [redBenMom1, redBenMom2], Aratio);
MRD = redBenMom * b * d^2  * fcd * 1000;