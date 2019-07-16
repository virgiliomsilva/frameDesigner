%% SHOT CALLER
clear
clc
fck = 30 ;
fyk = 400 ;
cover = .035 ;
% %% DC1
% %
% buildingName = 'regular_DC1' ;
% seismicCases = [24:31];
% folder = 'output\DC1' ;
% DC1frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
%% DC2
%
buildingName = 'regular_DC2' ;
seismicCases = [24:31];
folder = 'output\DC2' ;
DC2frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
% DC3

buildingName = 'regular_DC3' ;
seismicCases = [24:31];
folder = 'output\DC3' ;
DC3frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
%



%minimum stirrups to ensure proper bracing  - the method herein used is
%valid to equally spaced rebars
bracingDist = .3/2; %bracing distance on code
spaces = noRebar / 4;
inDist = (b - 2 * (cover + .02)); %distance from corner rebar to the middle of the section
switch noRebar
    case 8
        if inDist >= 2 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 3,:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3]),:);
        end
        
    case 12
        if inDist >= 3 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 4,:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4]),:);
        end
        
    case 16
        if inDist >= 4 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 5,:);
        elseif inDist < 2 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 3 4 5]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [3 4 5]),:);
        end
        
    case 20
        if inDist >= 5 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 6,:);
        elseif inDist < 2.5 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [2 4 6]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [4 6]),:);
        end
        
    case 24
        if inDist >= 6 * bracingDist
            shearReinforce = shearReinforce(shearReinforce(:,2) == 7,:);
        elseif inDist < 3 * bracingDist
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [3 4 5 6 7]),:);
        else
            shearReinforce = shearReinforce(ismember(shearReinforce(:,2), [4 5 6 7]),:);
        end
end