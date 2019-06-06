function [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = columnDesignDC3shear(fck, fyk , cover, Fz_Ed, b_input, h_input, longPhi)
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd; 

%%
sec_b = b_input;
h = h_input;
%%
z = h - 2 * (cover + .02);  %approximated 
Asw_s = max(sec_b * .08 * sqrt(fck) / fyk, Fz_Ed / (z * fywd * 1000 * 2.5)); %p.100 %9.4 '+' 9.5 EC2 %assuming cot(theta) = 2.5
for k = 1 : size(shearReinforce,1)
    if shearReinforce(k,4) - Asw_s > 0 %& (mod(shearReinforce(k,2),2) == mod(space + 1,2) | mod(shearReinforce(k,2) / 2,2) == mod(space + 1,2))
        diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
    else
        diffAuxS(k,1) = 1000;
    end
end
[minDiffS, minIndexS] = min(diffAuxS,[],1);


shearReinfPhi = shearReinforce(minIndexS, 1);
shearReinfLoops = shearReinforce(minIndexS, 2);
DC3_Spacing = min([h/4, 24 * shearReinfPhi /1000, 8 * longPhi]);

shearReinfSpac = min( shearReinforce(minIndexS, 3), DC3_Spacing);
shearReinfArea = shearReinforce(minIndexS, 4);