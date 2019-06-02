% function [sectionName, sectionType, sectionMaterials, sectionDimensions, reinforcementPatterns, ]= 
% sectionWriter(sectionName, fyk, fck, h, b, cover, ...
% rebarPhi, rebarQtd, stirPhi, stirBran, stirSpac)

clear
clc
tic

% columns = importdata('columnsEC2.mat');
% columnTypesAux = unique (columns(:,[2:9]),'rows');

sectionName =  'beam';
fyk = 500;
fck = 30;
h = .45;
b = .30;
cover = .035;
rebarPhi = 16
rebarQtd = 2
stirPhi = 8;
stirBran = 4;
stirSpac = .175;



% tf = strncmpi(s1,s2,4)

% if strncmpi(sectionName, 'column',3)
%     disp('accepted th column')
% else
%     disp('ohhh cryyy')
% end
%%

%% section material
switch fyk
    case {400, 500}
        steel = strcat('S', num2str(fyk));
    otherwise
        disp('Steel not supported');
end

switch fck
    case 25
        concrete = string ('C25/30');
    case 30
        concrete = string ('C30/37');
    otherwise
        disp('Concrete nod supported');
end

sectionMaterial = strcat(steel, "   ", concrete);
%% section dimensions
sectionDimensions = strcat(num2str(h), "   ", num2str(b), "   ", num2str(cover));
%% reinforcement patterns

%assuming beam

switch rebarPhi %& rebarQtd
    case {12, 16, 20, 25} %& {2, 3, 4}
        disp('FUCK YEAH!!')
    otherwise
        disp('YOU DUMB FUCK SON OF A BITCH!')
        
end

%% additional reinforcement

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
FRPstreng = string("wrapping(-)");
confinementFactor = string("auto135   (1   1)");
shearCapacity = string("yes   auto   Vc_Vp[+]   (1  1)");

toc