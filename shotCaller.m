%% SHOT CALLER
% data folder is the input folder; for each ductility class put the
% information on an appropriate folder and those are the names for
% buildingName#
clear 
clc
fck = 30 ;
fyk = 400 ;
cover = .035 ;
%% DC1
buildingName1 = 'regular_DC1' ;
nonSeismicCases = [22,23];
folder = ['output\' buildingName1];
mkdir(folder);
DC1frameDesigner(buildingName1, fck, fyk, cover, nonSeismicCases, folder);
%% DC2
buildingName2 = 'regular_DC2' ;
nonSeismicCases = [22,23];
seismicCases = [24:31];
folder = ['output\' buildingName2];
mkdir(folder);
DC2frameDesigner(buildingName2, fck, fyk, cover, seismicCases, nonSeismicCases, folder);
%% DC3
buildingName3 = 'regular_DC3' ;
nonSeismicCases = [22,23];
seismicCases = [24:31];
folder = ['output\' buildingName3];
mkdir(folder);
DC3frameDesigner(buildingName3, fck, fyk, cover, seismicCases, nonSeismicCases, folder, 'no');
% system('shutdown -s')