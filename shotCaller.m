%% SHOT CALLER
% data folder is the input folder; for each ductility class put the
% information on an appropriate folder and those are the names for
% buildingName#
[y,mo,d] = ymd(datetime);
[h,m,s] = hms(datetime);
start = {y,mo,d,h,m,s};
save start

clear; clc
fck = 30 ; %GPa
fyk = 400 ; %GPa
cover = .035 ; %m
%% DC1
% buildingName1 = 'tenerIT_DC1' ;
% seismicVerticalLoadCase = 40;
% nonSeismicCases = 23;%[22,23];
% folder = ['output\' buildingName1];
% mkdir(folder);
% DC1frameDesigner(buildingName1, fck, fyk, cover, seismicVerticalLoadCase,nonSeismicCases, folder);
%% DC2
buildingName2 = 'tenerIT_DC2' ;
seismicVerticalLoadCase = 40;
nonSeismicCases = 23;%[22,23];
seismicCases = [24:31];
folder = ['output\' buildingName2];
mkdir(folder);
DC2frameDesigner(buildingName2, fck, fyk, cover, seismicVerticalLoadCase, seismicCases, nonSeismicCases, folder);
%% DC3
buildingName3 = 'tenerIT_DC3' ;
seismicVerticalLoadCase = 40;
nonSeismicCases = 23;%[22,23];
seismicCases = [24:31];
folder = ['output\' buildingName3];
mkdir(folder);
DC3frameDesigner(buildingName3, fck, fyk, cover, seismicVerticalLoadCase, seismicCases, nonSeismicCases, folder, 'no');
%% DC3 w/ SLAB
buildingName4 = 'tenerIT_DC3slab' ;
seismicVerticalLoadCase = 40;
nonSeismicCases = 23;%[22,23];
seismicCases = [24:31];
slabTopReinf = 0.000314159; % top reinforcement slab m2/m       #10//.25
folder = ['output\' buildingName4];
mkdir(folder);
DC3frameDesignerWSlab (buildingName4, fck, fyk, cover, seismicVerticalLoadCase, seismicCases, nonSeismicCases, folder, 'no', slabTopReinf);
%%
[y,mo,d] = ymd(datetime);
[h,m,s] = hms(datetime);
finish = {y,mo,d,h,m,s};
save finish

system('shutdown -s')