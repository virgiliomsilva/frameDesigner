%% SHOT CALLER
clear
clc
fck = 30 ;
fyk = 400 ;
cover = .035 ;
%% DC1
%
buildingName = 'regular_DC1' ;
seismicCases = [24:31];
folder = 'output\DC1' ;
DC1frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
%% DC2
%
buildingName = 'regular_DC2' ;
seismicCases = [24:31];
folder = 'output\DC2' ;
DC2frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
%% DC3
%
buildingName = 'regular_DC3' ;
seismicCases = [24:31];
folder = 'output\DC3' ;
DC3frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder);
%%