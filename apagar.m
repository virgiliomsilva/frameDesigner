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

% a = minWidFind(89)
close(loading); loading = waitbar(0,'Initializing columns','Name', 'Step 3 of 6'); pause(1);
noStories = max(stories(:,1));
count = 0;
columns = [];
for i = 1 : size(barsOfColumns,1)
    % 1 design individually bars of a column
    barNames = []; % to append in the end!!!
    for j = 1 : noStories %design of all bars of a column!
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesign(:,1,1) == barName);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        for k = 1 : length(cases)
            N_axial = DataDesign(barIndex, 2, k);
            My_h = DataDesign(barIndex, 3, k);
            Mz_b = DataDesign(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1(k,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
        end
        [~, index] = max(mAux1(:,6));%best iteration for that bar
        bestIndi(j,:) = mAux1(index,:);%best individual bars of that column
    end
    
