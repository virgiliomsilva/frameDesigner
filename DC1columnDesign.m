%% DC1 COLUMNS DESIGN FUNCTION
% square columns with evenly distributed and spaced reinforcement along the
% four sides
function [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea,niter] = DC1columnDesign(fck, fyk , cover, N_Axial, My_h, Mz_b, givenWidth, givenLong,maxiter)%longRebarN, longRebarPh, 
%% APAGAR
% clear
% clc
% tic
% 
% fyk = 500;
% fck = 30;
% cover = .035;
% N_Axial = 900;
% My_h = 220;
% Mz_b = 180;
% givenWidth = 
%% INFO SELECTION
if (fyk == 500 & fck > 12 & fck < 50)
    abaco = importdata('info\abacus_REBAP83_C12_C50_S500.mat');
elseif (fyk == 400 & fck > 12 & fck < 50)
    abaco = importdata('info\abacus_REBAP83_C12_C50_S400.mat');
else
    error('Materials pair not supported!')
end

shearReinforce = importdata('info\steel_shear.csv');
if exist('givenLong', 'var')
%     if size(givenLong,1) > 1
    longReinforce = givenLong;
%     else
%         longReinforce = columnComp(givenLong);
%     end
else
    longReinforce  = importdata('info\steel_columnEC8.csv'); 
    longReinforce = longReinforce(:,[1:3]);
end

fcd = fck / 1.5; 
fyd = fyk / 1.15;
dMax = .03; %maximum aggregate dimension
%% DIMENSIONS & LONG REBAR
if exist ('givenWidth','var')
    h = givenWidth;
else
    h = .2; 
end
b = h;

incr = .05;

diffe = -1; areaRebar = 2; AsMax = 1; %values to make while loop start
niter=0;

if ~exist('maxiter','var')
    maxiter=10;
end
    
while (diffe < 0 | areaRebar > AsMax %| AsMin > areaRebar)%& niter < maxiter
    redAxial   = N_Axial / (b * h * fcd * 1000);
    redBenMom1 = Mz_b / (b * h^2 * fcd * 1000);
    redBenMom2 = My_h / (h * b^2 * fcd * 1000);
    bigRedBenMom = max(redBenMom1, redBenMom2);
    smallRedBenMom = min(redBenMom1, redBenMom2);
    redBenMomRatio = bigRedBenMom / smallRedBenMom;
    
    reinfPerc = abaco(redBenMomRatio, redAxial, smallRedBenMom);
    AsAbacus = reinfPerc * b * h * fcd / fyd;
    
    AsMin = max([.1 * N_Axial / (fyd * 1000), .01 * b * h]);% first condition?? second eurocode 8 otherwise 0.2%
    reinfArea = max([AsMin, AsAbacus]);
        
    for j = 1 : size(longReinforce,1)
        if longReinforce(j,3) - reinfArea > 0
            diffAux(j) = longReinforce(j,3) - reinfArea;
        else
            diffAux(j) = 1000;
        end
    end

    [~, minIndex] = min(diffAux);
    noRebar = longReinforce(minIndex,2);
    phiRebar = longReinforce(minIndex,1);
    areaRebar = longReinforce(minIndex,3);
    
    diffe = b - 2 * (cover + .01) - (phiRebar/1000 * (noRebar/4 + 1)) - (max([.02, phiRebar/1000, dMax+.005]) * (noRebar/4));
    AsMax = .04 * h * b; %EC2 & EC8
    %AsMin = .01 * h * b; %EC8
    
    h = h + incr; b = h; niter = niter+1;
    if niter == maxiter
        error('Max iter reached')
    end
end
sec_h = h - incr;    sec_b = sec_h ;

% M_Rd 
reinfPercFin = areaRebar * (b*h);%fyd / (b * h * fcd);
redAxialFin = N_Axial / (b * h * fcd * 1000);
diff = [];
for i = 0 : .005 : .5
    reinfPerc = abaco(1, redAxialFin, i);
    diff = [diff; [i , abs(reinfPercFin - reinfPerc)]];
end

[~, index1] = min(diff(:, 2));
redBenMom = (index1 - 1) * .005;
M_Rd = redBenMom * h * b^2 * fcd * 1000 ; 

%N_rd
alpha = 6;
N_rd = fcd * (sec_h * sec_b + areaRebar * alpha);

%% STIRRUPS
%maximum spacing
max_spacing = min([15 * phiRebar/1000, b, .3]) * .6; %evaluating the shear spacing on the critical sections
shearReinforce(shearReinforce(:,3) > max_spacing, :) = [];

%possible loops due to rebar configuration - deprecated due to overlapping
%to code block "minimum stirrups to ensure proper bracing"
% switch noRebar
%     case 8
%         shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3]),:)
%     case 12
%         shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4]),:)
%     case 16
%         shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3 4 5]),:)
%     case 20
%         shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4 6]),:)
%     case 24
% end

%minimum stirrups to ensure proper bracing  - the method herein used is
%valid to equally spaced rebars
bracingDist = .15; %bracing distance on code
spaces = noRebar / 4;
halfDist = (b - 2 * (cover + .02)) / 2; %distance from corner rebar to the middle of the section
switch noRebar
    case 8
        if halfDist >= bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 3,:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3]),:);
        end
        
    case 12
        if halfDist >= 1.5 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 4,:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4]),:);
        end
        
    case 16
        if halfDist >= 2 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 5,:);
        elseif halfDist < 1 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3 4 5]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [3 4 5]),:);
        end
        
    case 20
        if halfDist >= 2.5 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 6,:);
        elseif halfDist < 1.25 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4 6]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [4 6]),:);
        end  
            
    case 24 
        if halfDist >= 3 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 7,:);
        elseif halfDist < 1.5 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [3 4 5 6 7]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [4 5 6 7]),:);
        end
end

[minDiffS, minIndexS] = min(shearReinforce(:,4));
shearReinfPhi = shearReinforce(minIndexS, 1);
shearReinfSpac = shearReinforce(minIndexS, 3);
shearReinfLoops = shearReinforce(minIndexS, 2);
shearReinfArea = shearReinforce(minIndexS, 4);