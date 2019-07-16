function [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = DC3beamDesignMidShear(fck, fyk , cover, Fz_Ed, sec_b, sec_h, longReinfNo)
%% INFO SELECTION

shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5;
fyd = fyk / 1.15;
fywd = fyd;

%% STIRRUPS

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

redVed = Fz_Ed / (sec_b * z * 1000);
redVrd25 = (.6/4.35) * fck * (1 - fck/250);
redVrd1 = (.6/3) * fck * (1 - fck/250);

if redVrd25 > redVed
    Asw_s25 = redVed * sec_h / (fywd * 2.5);
    Asw_sMin = .08*sqrt(fck)/fyk  * sec_b;
    Asw_s = max([Asw_sMin, Asw_s25]);
    theta = acot(2.5);
elseif redVrd1 > redVed
    theta = .5 * asin(redVed / (.20 * fck * (1 - fck/250)));
    Asw_stheta = redVed * sec_h / (fywd * theta);
    Asw_sMin = .08*sqrt(fck)/fyk  * sec_b;
    Asw_s = max([Asw_sMin, Asw_stheta]);
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

V_Rd = min(shearReinfArea * z * (fywd * 1000) * cot(theta) / shearReinfSpac, 1 * sec_b * z * .6 *(1 - fck/250) * fcd  * 1000/ (cot(theta) + tan(theta)));
end