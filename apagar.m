% clear
% clc
% N_axial = 804;% max(DataDesign(DataDesign(:,1,1) == toGive(j,1), 2, seismicCasesIdx));
% My_h = 283;
% Mz_b = 0;%My_h / 3;
% % [columnsRow,~] = find(columns(:,1) == toGive(j,1))
% gWidth = .35;%columns(columnsRow, 2);
% % input1 = columns(columnsRow, 4);
% % input2 = columns(columnsRow, 5);
% % g_reinf = columnComp(input1, input2, 'EC8', 'yes');
% fck=30; fcd= 20;
% fyk=400;fyd= fyk/1.15;
% cover= .035;
% g_reinf = [16,12,0.00241274300000000;16,16,0.00321699100000000;20,12,0.00376991100000000;20,16,0.00502654800000000];
% 
% 
% [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, g_reinf);% [20,16,0.00502654800000000]);%, gWidth);
% % b = .3; h = b;
% % abaco = importdata('info\abacus_REBAP83_C12_C50_S400.mat');
% % redAxial   = N_axial / (b * h * fcd * 1000);
% % redBenMom1 = Mz_b / (b * h^2 * fcd * 1000);
% % redBenMom2 = My_h / (h * b^2 * fcd * 1000);
% % bigRedBenMom = max(redBenMom1, redBenMom2);
% % smallRedBenMom = min(redBenMom1, redBenMom2);
% % redBenMomRatio = smallRedBenMom / bigRedBenMom;
% % % redBenMomRatio
% % % redAxial
% % % smallRedBenMom
% % reinfPerc = abaco(redBenMomRatio, redAxial, bigRedBenMom);
% % function [res] = apagar (a,b)
% res = a * b ;

[newComp] = columnComp(16, 16, 'EC8', 'yes')
