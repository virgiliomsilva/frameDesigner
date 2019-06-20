%% DC1frameDesigner
% function [] = DC1frameDesigner()

clear
clc
tic
disp('Started');
buildingName = 'regular';
fnData = ['data\' buildingName '\datasetcalcada.csv'] ;
fnNodes = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;
fck = 30;
fyk = 400;
cover = .035;

% [barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories, nodes, cases] = ...
%     dataTransformer (fnData, fnElement, fnNodes);
[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesign, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);
%%
beams = [];
loading = waitbar(0,'Initializing beams'); pause(.5);
% for i = 1 : length(beamDesiOrd)
%     barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
%     
%     mAux = [];
%     for j = 1 : length(cases)
%         M_Ed = DataDesign(barIndex, 5, j);
%         [sec_h, sec_b, longReinfNo, longReinfPhi, ~, roMinCondition, M_R] = DC1beamDesign(fck, fyk , cover, M_Ed, 0);
%         mAux(j,:) = [sec_h, sec_b, longReinfNo, longReinfPhi, roMinCondition, M_R];
%     end
%     [M_rd, conIndex] = max(mAux(:,6)); %best M_rd 
% 
%     sAux = [];
%     for j = 1 : length(cases)
%         Fz_Ed = DataDesign(barIndex, 4, j);
%         given_h = mAux(conIndex, 1);
%         given_b = mAux(conIndex, 2);
%         longRebarN = mAux(conIndex, 3);
%         longRebarPh = mAux(conIndex, 4);
%         [~,~,~,~,~,~,~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it] = DC1beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
%         sAux(j,:) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it];
%     end
%     [V_Rd, conIndex2] = max(sAux(:,4)); %best F_rd
%     
%     shearPhi = sAux(conIndex2, 1);
%     shearSpac = sAux(conIndex2, 2);
%     shearLegs = sAux(conIndex2, 3);
% 
%     beams(end+1,:) = [DataDesign(i,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, shearPhi, shearSpac, shearLegs, V_Rd];
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress'); %pause(.5);
% end

beams = importdata('beams.mat');
disp(['Time to design beams: ' num2str(toc)])
%% tic
close(loading); loading = waitbar(0,'Initializing columns'); pause(1);
noStories = max(stories(:,1));
count = 0;
countf =0;
% columns = [];
% for i = 1 : size(barsOfColumns,1) 
%     % 1 design individually bars of a column
%     barNames = []; % to append in the end!!!
%     for j = 1 : noStories %design of all bars of a column!
%         barName = barsOfColumns(i,j); barNames = [barNames; barName];
%         barIndex = find(DataDesign(:,1,1) == barName);
%         for k = 1 : length(cases)
%             N_axial = DataDesign(barIndex, 2, k);
%             My_h = DataDesign(barIndex, 3, k);
%             Mz_b = DataDesign(barIndex, 4, k);
%             [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b);
%             mAux1(k,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
%         end
%         [~, index] = max(mAux1(:,6));%best iteration for that bar
%         bestIndi(j,:) = mAux1(index,:);%best individual bars of that column
%     end
% %
% % the snippet
% %     
%     % 3 design from bottom to top
%     %bottom bar
%     bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration    %gWidth = floor(bigOrigWidth * .9 * 20)/20;
%     gWidth = bigOrigWidth;
%     barName = barsOfColumns(i,1);
%     barIndex = find(DataDesign(:,1,1) == barName);
%     for p = 1 : length(cases)
%         N_axial = DataDesign(barIndex, 2, p);
%         My_h = DataDesign(barIndex, 3, p);
%         Mz_b = DataDesign(barIndex, 4, p);
%         [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth);
%         ratio = bestIndi(1,6) / reinfPercFin;
%         mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];%mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, ratio, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
%     end
%     [~, index] = max(mAux3(:,6));%best iteration for that bar
%     auxColumns = mAux3(index,:);
% %     bestIndi2(1,:,k+1) = mAux3(index,:);
%             
%     %other bars
%     gWidth = floor(mAux3(index,1) * .9 * 20)/20;    %givenLong = mAux3(index, 5);
%     for j = 2 : noStories
%         barName = barsOfColumns(i,j);
%         barIndex = find(DataDesign(:,1,1) == barName);
%         for p = 1 : length(cases)
%             N_axial = DataDesign(barIndex, 2, p);
%             My_h = DataDesign(barIndex, 3, p);
%             Mz_b = DataDesign(barIndex, 4, p);
%             givenLong = columnComp(mAux3(index, 5),'EC8');
%             [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
%             %ratio = bestIndi(j,6) / reinfPercFin;
%             mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];%mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, ratio, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
%         end
%         [~, index] = max(mAux3(:,6));%best iteration for that bar
%         auxColumns(j,:) = mAux3(index,:);
% %         bestIndi2(j,:,k+1) = mAux3(index,:);%best individual bars of that column with a seed provided by the parent bar
%         try
%         gWidth = floor(auxColumns(j+1,1) * .9 * 20)/20;
%         givenLong = mAux3(index, 6);
%         end
%     end
%     
% %     possDesigns = [];
% %     for z = 1 : (noStories + 1) % check if superior dimesions are the same or smaller, if not, kick it
% %         agr = issorted(bestIndi2(:,1,z),'descend');
% %         if agr == 1
% %             possDesigns = [possDesigns, z];
% %         end
% %     end
% %     validDesigns = bestIndi2(:,:,possDesigns);
%     %check if dimensions agree
%     sorted = issorted(auxColumns(:,1),'descend');
%     if sorted == 0
%         [row, col] = size(auxColumns)
%         auxColumns = zeros(row, col);
%         error('Columns with strange dimensions!!')
%         %TODO add here the other dimensioning method, if that ~fails.. well
%     end
%     
%     columns = [columns; [barNames, auxColumns]];
%     
%     count = count + 1;
% %     disp('');
% %     disp(['Finished the ' num2str(count) ' column out of ' num2str(size(barsOfColumns,1))]);
% %     disp(['Calls of beamDesign function: ' num2str(countf)]);
%     toc
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress');% pause(1);
%    % if count == 2; break; end
% end
columns = importdata('columns.mat');
%% 1.3 comparison
close(loading); loading = waitbar(0,'Initializing bending comparisons'); pause(1);
increNeed = [];
for i = 1 : size(nodes,1)
    if nodes(i, 5) == 0 | nodes(i, 5) == noStories
        continue %skip base points and on the top floor the verification is not needed
    end
    
    %get bars on that point
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 1); barsX = element(row,1);
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 2); barsY = element(row,1);
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 3); barsZ = element(row,1);
    
    %get design bending moments on those bars
    bendRdX = 0;
    for j = 1 : length(barsX)
        [beamRow, ~] = find(beams(:,1) == barsX(j));
        bendRdX = bendRdX + beams(beamRow, 6);
    end

    bendRdY = 0;
    for j = 1 : length(barsY)
        [beamRow, ~] = find(beams(:,1) == barsY(j));
        bendRdY = bendRdY + beams(beamRow, 6);
    end
    
    bendRdZ = 0;
    for j = 1 : length(barsZ)
        [columnRow, ~] = find(columns(:,1) == barsZ(j));
        bendRdZ = bendRdZ + columns(beamRow, 6);
    end
    
    %do the sums and comparisons
    ratioXZ = bendRdX / bendRdZ;
    if ratioXZ < 1.3
        auxBendMomX = bendRdZ * 1.3;
    else
        auxBendMomX = 0;
    end
    
    ratioYZ = bendRdY / bendRdZ;
    if ratioYZ < 1.3
        auxBendMomY = bendRdY * 1.3;
    else
        auxBendMomY = 0;
    end
    
    increNeed(end+1,:) = [nodes(i,1), max(auxBendMomX, auxBendMomY)];
    clear auxBendMomX auxBendMomY
    loading = waitbar(i / (size(nodes,1) * (noStories - 1)),loading,'Comparisons progress');
end

% organise the values
% get the bars of that point
% get both reinforcement patterns possible compatiblities
% find the mutual 
% for that size get the needed bending 
% check if dimensions have changed


%for each column find the significant value to increase from top to
%bottom!!




%toSeismo(columns, beams, nodes, element, stories, fck, fyk, cover)
close(loading);
disp('Finished');
toc