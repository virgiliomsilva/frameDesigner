%% hardest way to get the MRD of a beam!!
function [MRD] = MrdBeam (fck, fyk, cover, b, h, A1, A2) 
if (fyk == 500 & fck > 12 & fck < 50)
    abaco = importdata('info\abacus_C12_50A500_allBend.mat');
elseif (fyk == 400 & fck > 12 & fck < 50)
    abaco = importdata('info\abacus_C12_50A400_allBend.mat');
else
    error('Materials pair not supported!')
end

fyd = fyk / 1.15;
fcd = fck / 1.5;

Aratio = min([A1, A2]) / max([A1, A2]);

d = h - (cover + .02);
reinfPerc = max([A1, A2]) * fyd / (b * d * fcd);

Ind = find(abaco(1,:) > Aratio, 1); 
% if
    
for i = 2 : size(abaco,1)
     auxR(i) = ((abaco(i, Ind) - abaco(i, Ind-1)) / (abaco(1, Ind) - abaco(1, Ind-1))) * (Aratio - (abaco(1, Ind-1))) + 
end