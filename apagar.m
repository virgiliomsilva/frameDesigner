% clear
% clc
% tic
% table = importdata('info\steel_ColumnEC8.csv');
% table = table(:, [1:3]);
% 
% fck = 30;
% fyk = 400;
% cover = .035;
% N_axial = 0;
% My_h = 0;
% Mz_b = 0;
% gWidth = .3;
% 
% 
% res = [];
% loading = waitbar(0,'Initializing'); pause(1);
% for i = 1 : size(table,1)
%     givenLong = table(i,:);
%     
%     [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
% 
%     res = [res; [sec_h, sec_b, noRebar, phiRebar, reinfPercFin, M_Rd]];
%     
%     loading = waitbar( i / size(table,1),loading,'Progress');
% end
% 
% res = sortrows(res,6);
% pause(2);
% close(loading);
% toc
table = [1 2 3; 1 3 2]
[row, col] = find(table == 3)