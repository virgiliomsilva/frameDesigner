function [result]= sectionWriter(sectionName, fyk, fck, h, b, cover, ...
rebarPhi, rebarQtd, stirPhi, stirBran, stirSpac)
% sectionName =  'beamntttt';
% fyk = 400;
% fck = 30;
% h = .45;
% b = .30;
% cover = .035;
% rebarPhi = 22.5;
% rebarQtd = 8;
% stirPhi = 8;
% stirBran = 4;
% stirSpac = .175;
%% section name
originalSectionName = sectionName;

if strncmpi(sectionName, 'column',4)
    sectionName = 'column' ;
elseif strncmpi(sectionName, 'beam',4)
    sectionName = 'beam' ;
end
%% section type
sectionType = 'rcrs';
%% section material
switch fyk
    case {400, 500}
        steel = strcat('S', num2str(fyk));
    otherwise
        steel = "Steel not supported";
end

switch fck
    case 25
        concrete = string ('C25/30');
    case 30
        concrete = string ('C30/37');
    otherwise
        concrete = "Concrete not supported";
end

sectionMaterial = strcat(steel, "   ", concrete);
%% section dimensions
sectionDimensions = strcat(num2str(h), "   ", num2str(b), "   ", num2str(cover));
%% reinforcement patterns & additional reinforcement
switch sectionName
    case 'beam'
        switch rebarPhi
            case {12, 16, 20, 25}
                switch rebarQtd
                    case 2
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(0@16mm)   left_right_sides(0@16mm)');
                        addReinf = {};
                    case 3
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(2@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = {};
                    case 4
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(2@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = {};
                    case 5
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(2@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, rebarPhi, rebarQtd, stirPhi);
                    case 6
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(2@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, rebarPhi, rebarQtd, stirPhi);
                    case 8
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(4@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, rebarPhi, rebarQtd, stirPhi);
                    case 9
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(6@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, rebarPhi, rebarQtd, stirPhi);
                    case 10
                        reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(6@', num2str(rebarPhi),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, rebarPhi, rebarQtd, stirPhi);
                end

            case {14, 18 22.5}
                newRebarPhiOut = ceil(rebarPhi + 2);
                newRebarPhiIn = floor(rebarPhi - 2);
                switch rebarQtd
                    case 4
                        reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(4@', num2str(newRebarPhiIn),'mm)   left_right_sides(0@16mm)');
                        addReinf = {};
                    case 6
                        reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(2@', num2str(newRebarPhiOut),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, newRebarPhiIn, rebarQtd, stirPhi);
                    case 8
                        reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(4@', num2str(newRebarPhiOut),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, newRebarPhiIn, rebarQtd, stirPhi);
                    case 10
                        reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(6@', num2str(newRebarPhiOut),'mm)   left_right_sides(0@16mm)');
                        addReinf = addiReinfBeam (b, h, cover, newRebarPhiIn, rebarQtd, stirPhi);
                end


            case {14.7, 18.7, 23.3}
                newRebarPhiOut = ceil(rebarPhi + 1);
                newRebarPhiIn = floor(rebarPhi - 2.5);
                reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(2@', num2str(newRebarPhiIn),'mm)   left_right_sides(0@16mm)');
                addReinf = {};

            case {18.4, 23}
                newRebarPhiOut = ceil(rebarPhi + 1.5);
                newRebarPhiIn = floor(rebarPhi - 2.2);
                reinfPatt = strcat('corners(4@', num2str(newRebarPhiOut), 'mm)   top_bottom_sides(2@', num2str(newRebarPhiOut),'mm)   left_right_sides(0@16mm)');
                addReinf = addiReinfBeam (b, h, cover, newRebarPhiIn, rebarQtd, stirPhi);

            otherwise
                reinfPatt = "Reinforcement pattern not supported";
                addReinf = reinfPatt;
        end

    case 'column'
        switch rebarPhi
            case {12, 16, 20, 25, 32}
                noRebarSide = (rebarQtd - 4) / 2; 
                reinfPatt = strcat('corners(4@', num2str(rebarPhi), 'mm)   top_bottom_sides(', num2str(noRebarSide),'@', num2str(rebarPhi),'mm)   left_right_sides(', num2str(noRebarSide),'@', num2str(rebarPhi),'mm)');
            case 13.1
                reinfPatt = 'corners(4@16mm)   top_bottom_sides(6@12mm)   left_right_sides(6@12mm)';
            case 13.5
                reinfPatt = 'corners(4@16mm)   top_bottom_sides(4@12mm)   left_right_sides(4@12mm)';
            case 14.1
                reinfPatt = 'corners(4@16mm)   top_bottom_sides(2@12mm)   left_right_sides(2@12mm)';
            case 17.1
                reinfPatt = 'corners(4@20mm)   top_bottom_sides(6@16mm)   left_right_sides(6@16mm)';
            case 17.4
                reinfPatt = 'corners(4@20mm)   top_bottom_sides(4@16mm)   left_right_sides(4@16mm)';
            case 18.1
                reinfPatt = 'corners(4@20mm)   top_bottom_sides(2@16mm)   left_right_sides(2@16mm)';
            case 21.4
                reinfPatt = 'corners(4@25mm)   top_bottom_sides(6@20mm)   left_right_sides(6@20mm)';
            case 21.8
                reinfPatt = 'corners(4@25mm)   top_bottom_sides(4@20mm)   left_right_sides(4@20mm)';
            case 22.6
                reinfPatt = 'corners(4@25mm)   top_bottom_sides(2@20mm)   left_right_sides(2@20mm)';
            case 26.9
                reinfPatt = 'corners(4@32mm)   top_bottom_sides(6@25mm)   left_right_sides(6@25mm)';
            case 27.5
                reinfPatt = 'corners(4@32mm)   top_bottom_sides(4@25mm)   left_right_sides(4@25mm)';
            case 28.7
                reinfPatt = 'corners(4@32mm)   top_bottom_sides(2@25mm)   left_right_sides(2@25mm)';
        end
        addReinf = {};
end
%% transverse reinforcement
switch sectionName
    case 'column'
        reinfHori = stirBran;
        reinfVert = stirBran;
    case 'beam'
        reinfHori = 2;
        reinfVert = stirBran;
end

transverseReinforcement = strcat("(", num2str(reinfHori), "-", num2str(reinfVert), ")", "@", num2str(stirPhi), "mm/", num2str(stirSpac));
%% FRP strengthening & confinement factors & shear capacity
FRPstreng = "wrapping(-)" ;
confinementFactor = "auto135   (1   1)" ;
shearCapacity = "yes   auto   Vc_Vp[+]   (1  1)" ;
%% 
result = {originalSectionName, sectionType, sectionMaterial, sectionDimensions, ...
    reinfPatt, addReinf, transverseReinforcement, FRPstreng, confinementFactor, shearCapacity};