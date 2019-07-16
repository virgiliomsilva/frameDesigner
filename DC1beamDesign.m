%% DC1 BEAM DESIGN FUNCTION
% beams with a maximum dimension of width over height of 80% with equal
% reinforcement on top and bottom it can design only stirrups given the
% longitudinal rebar
function [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition] = DC1beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longReinfN, longReinfPh)
sCondition = 0;
%% INFO SELECTION
if (fyk == 500 & fck > 12 & fck < 50)
    abaco = importdata('info\abacusC12_50S500A1.csv');
elseif (fyk == 400 & fck > 12 & fck < 50)
    abaco = importdata('info\abacusC12_50S400A1.csv');
else
    error('Materials pair not supported!')
end

longReinforce  = importdata('info\steel_beam.csv'); longReinforce = longReinforce(:, [1:3]);
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5;
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd;
dMax = .03;
%% DIMENSIONS & LONG REBAR
h = .25;
b = .2;
incr = .05;
count = 0; M_Rd = 0;
while M_Rd < M_Ed || bOh > .8 || longReinfArea > AsMax
    
    %long rebar starting values
    d = h - (cover + .02);  %approximated
    redBendMom = M_Ed / (b * d^2 * fcd * 1000);
    reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]),'linear', 1);
    reinfAreaAux = reinfPerc * b * d * fcd / fyd;
    [reinfArea, roMinCondition] = max([reinfAreaAux, .26*fctm*b*d/fyk, .0013*b*d, .5*fctm/fyk*b*d]);
    % condition 1 - given by the abacus
    % condition 2 - EC2 9.1
    % condition 3 - EC2 9.2.1.1 (9.1)
    % condition 4 - new EC8 10.5.4.2 (10.2)
    
    %long rebar iterations
    if reinfPerc > 0 && reinfPerc < .8
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfArea > 0
                diffAuxL(j) = longReinforce(j,3) - reinfArea;
            else
                diffAuxL(j) = 1;
            end
            [minDiff, minIndex] = min(diffAuxL);
            longReinfNo = longReinforce(minIndex,2); %for each zone (comp or tension)
            longReinfPhi = longReinforce(minIndex,1);
            longReinfArea = longReinforce(minIndex,3);
        end
        
        if longReinfNo <= 4
            spaces = longReinfNo - 1;
            clearance = spaces * max([longReinfPhi/1000, dMax + .005, .02]);
        else
            spaces = ceil(longReinfNo / 2) - 1;
            clearance = spaces * max([longReinfPhi/1000, dMax + .005, .02]);
        end
        
        diff = b - 2 * (cover + .01) - clearance - (spaces + 1) * (longReinfPhi / 1000);
        while diff < 0
            b = b + incr;
            diff = b - 2 * (cover + .01) - clearance - (spaces + 1) * (longReinfPhi / 1000);
        end
        
        d = h - (cover + .02);
        reinfPerc = longReinfArea * fyd / (b * d * fcd);
        redBendMom = interp1(abaco(:,2).', abaco(:,1).', max([.004, reinfPerc]), 'linear', 0);
        M_Rd = redBendMom * b * d^2 * fcd * 1000 ;
    end
    
    sec_h = h ;
    sec_b = b ;
    AsMax = .04 * b * h;
    bOh = b / h ;
    % values for the next iterarion on dimensions
    h = .25 + (count + 1) * incr;
    b = max([floor(.3 * h * 20) / 20, .2]);
    count = count+1;
end
%% STIRRUPS
if exist('given_b', 'var'); sec_b = given_b; end
if exist('given_h', 'var'); sec_h = given_h; end
if exist('longReinfN', 'var'); longReinfNo = longReinfN; end
if exist('longReinfPh', 'var'); longReinfPhi = longReinfPh; end

%possible loops due to rebar configuration
shearReinforce(ismember(shearReinforce(:,2), [6 7]), :) = [] ;%no beam configuration support 6 or 7 stirrups
switch longReinfNo
    case 2
        shearReinforce(shearReinforce(:,2) ~= 2, :) = [];
    case {3, 5, 6, 9, 10}
        shearReinforce(shearReinforce(:,2) == 4, :) = [];
        shearReinforce(shearReinforce(:,2) == 5, :) = [];
    case {4, 8}
        shearReinforce(shearReinforce(:,2) == 3, :) = [];
        shearReinforce(shearReinforce(:,2) == 5, :) = [];
end

%maximum spacing longitudinal
d = sec_h - (cover + .02);  %approximated
max_spacing = .75 * d; %sl, max
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
else
    display('-----------------------------------------------------------');
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