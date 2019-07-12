% 
% clear
% clc
% fck = 30 ;
% fyk = 400;
% cover = .035 ;
% % dataTransformer está condicionado!! apagar linha 48
% %% DC1
% %
% % buildingName = 'regular_DC1' ;
% % seismicCases = [24:31];
% % DC1frameDesigner(buildingName, fck, fyk, cover, seismicCases);
% %% DC2
% %
% % buildingName = 'regular_DC2' ;
% % seismicCases = [24:31];
% % DC2frameDesigner(buildingName, fck, fyk, cover, seismicCases);
% %% DC3
% %
% % buildingName = 'regular_DC3' ;
% % seismicCases = [24:31];
% % DC3frameDesigner(buildingName, fck, fyk, cover, seismicCases);
% %% EC2
% %
fck = 30;
fyk = 400;
cover = .035;
N_Axial = 10;
My_h = 0;
Mz_b = 0;
givenWidth = .25;
givenLong = [18.1,8,0.00206088500000000]; 
%[18.1,8,0.00206088500000000]
% [16,12,0.00241274300000000]
[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_Axial, My_h, Mz_b, givenWidth);%, givenLong);%, Vshear, given_h, longReinfN, longReinfPh)