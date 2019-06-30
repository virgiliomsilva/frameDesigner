%% SHOT CALLER
clear
clc
fck = 30 ;
fyk = 400;
cover = .035 ;
% dataTransformer está condicionado!! apagar linha 48
%% DC1
%
buildingName = 'regular_DC1' ;
seismicCases = [24:31];
DC1frameDesigner(buildingName, fck, fyk, cover, seismicCases);
%% DC2
%
buildingName = 'regular_DC2' ;
seismicCases = [24:31];
DC2frameDesigner(buildingName, fck, fyk, cover, seismicCases);
%% DC3
%

%% EC2
%
