%% DC2 COLUMNS DESIGN FUNCTION
% square columns with evenly distributed and spaced reinforcement along the
% four sides
function [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC2columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar)
    %% INFO SELECTION
    if ~((fyk == 500 && fck > 12 && fck < 50) || (fyk == 400 && fck > 12 && fck < 50))
        error('Materials pair not supported!')
    end

    shearReinforce = importdata('info\steel_shear.csv');

    %% STIRRUPS
    %maximum spacing
    max_spacing = min([15 * phiRebar/1000, sec_h, .3]);
    shearReinforce(shearReinforce(:,3) > max_spacing, :) = [];

    %minimum stirrups to ensure proper bracing  - the method herein used is
    %valid to equally spaced rebars
    bracingDist = .2/2; %bracing distance on code EC8
    halfDist = (sec_h - 2 * (cover + .02)) / 2; %distance from corner rebar to the middle of the section
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

    [~, minIndexS] = min(shearReinforce(:,4));

    shearReinfPhi = shearReinforce(minIndexS, 1);
    shearReinfSpac = shearReinforce(minIndexS, 3);
    shearReinfLoops = shearReinforce(minIndexS, 2);
    shearReinfArea = shearReinforce(minIndexS, 4);
end