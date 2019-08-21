function [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = DC2beamDesignMidShear(fck, fyk , cover, Fz_Ed, sec_b, sec_h, longReinfNo, longReinfPhi)
%% INFO SELECTION
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5;
fyd = fyk / 1.15;
fywd = fyd;

Fz_Ed = abs(Fz_Ed);
%% STIRRUPS
longReinfArea = longReinfNo * pi * (longReinfPhi/2000)^2;

%possible loops due to rebar configuration
shearReinforce(ismember(shearReinforce(:,2), [6 7]), :) = [] ;%no beam configuration support 6 or 7 stirrups
switch longReinfNo
    case 2
        shearReinforce(shearReinforce(:,2) ~= 2, :) = [];
    case {3, 5, 6}
        shearReinforce(shearReinforce(:,2) == 4, :) = [];
        shearReinforce(shearReinforce(:,2) == 5, :) = [];
    case {4, 8}
        shearReinforce(shearReinforce(:,2) == 3, :) = [];
        shearReinforce(shearReinforce(:,2) == 5, :) = [];
end

%maximum spacing longitudinal
d = sec_h - (cover + .02);  %approximated
max_spacing = min([.75 * d, .3]);
shearReinforce(shearReinforce(:,3) > max_spacing, :) = [];

%maximum spacing between legs
st = min(.75 * d, .6); %p.176
extLegsDist = sec_b - 2 * (cover + .005);
minNoLegs = floor(extLegsDist / st) + 2;
shearReinforce(shearReinforce(:,2) < minNoLegs, :) = [];

% solution calculation
z = sec_h - 2 * (cover + .02);  %approximated

Asw_sMin = .08*sqrt(fck)/fyk  * sec_b;
Vrds = Asw_sMin * z * fywd * 1000 * 2.5;

kFactor = min(1 + sqrt(.2/d), 2);
roL = min(longReinfArea / (sec_b * sec_h), .02);
Vrdc1 = .12 * kFactor * (100 * roL * fck)^(1/3) ;
Vrdc2 = .035 * kFactor ^ 1.5 * sqrt(fck);
Vrdc = max(Vrdc1, Vrdc2) * sec_b * sec_h;

VRd = max(Vrdc, Vrds);

cotTheta = 2.5;
miu = .6 * (1 - fck/250);
if VRd > Fz_Ed
    Asw_s = Asw_sMin;
else
    VRd_max = sec_b * z * miu * fcd * 1000 / (2.5 + 1/2.5); 
    if VRd_max > Fz_Ed
        Asw_s = Fz_Ed / (z * fywd * 1000 * 2.5);
    else
        sol = solve(sec_b * z * miu * fcd * 1000 / (x + 1/x) == Fz_Ed, x);
        cotTheta = symsToCotTheta(sol);
        Asw_s = Fz_Ed / (z * fywd * 1000 * cotTheta);
    end
end

%minimum optimization for steel usage (p.82)
diffAuxS = zeros(size(shearReinforce,1),1);
for k = 1 : size(shearReinforce,1)
    if shearReinforce(k,4) - Asw_s > 0
        diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
    else
        diffAuxS(k,1) = 1000;
    end
end

%output
[~, minIndexS] = min(diffAuxS,[],1);
shearReinfPhi = shearReinforce(minIndexS, 1);
shearReinfSpac = shearReinforce(minIndexS, 3);
shearReinfLoops = shearReinforce(minIndexS, 2);
shearReinfArea = shearReinforce(minIndexS, 4);

VRd_max = sec_b * z * miu * fcd * 1000 / (cotTheta + 1/cotTheta);
Vrds = shearReinfArea * z * fywd * 1000 * cotTheta;
V_Rd = min (VRd_max, Vrds);end