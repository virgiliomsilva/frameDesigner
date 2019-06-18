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
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
    
    mAux = [];
    for j = 1 : length(cases)
        M_Ed = DataDesign(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, ~, roMinCondition, M_R] = DC1beamDesign(fck, fyk , cover, M_Ed, 0);
        mAux(j,:) = [sec_h, sec_b, longReinfNo, longReinfPhi, roMinCondition, M_R];
    end
    [M_rd, conIndex] = max(mAux(:,6)); %best M_rd 

    sAux = [];
    for j = 1 : length(cases)
        Fz_Ed = DataDesign(barIndex, 4, j);
        given_h = mAux(conIndex, 1);
        given_b = mAux(conIndex, 2);
        longRebarN = mAux(conIndex, 3);
        longRebarPh = mAux(conIndex, 4);
        [~,~,~,~,~,~,~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it] = DC1beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
        sAux(j,:) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it];
    end
    [V_Rd, conIndex2] = max(sAux(:,4)); %best F_rd
    
    shearPhi = sAux(conIndex2, 1);
    shearSpac = sAux(conIndex2, 2);
    shearLegs = sAux(conIndex2, 3);

    beams(end+1,:) = [DataDesign(i,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, shearPhi, shearSpac, shearLegs, V_Rd];
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress'); %pause(.5);
end

disp(['Time to design beams: ' num2str(toc)])
% tic
loading = waitbar(0,loading,'Initializing columns'); pause(1);
noStories = max(stories(:,1));
count = 0;
countf =0;
columns = [];
for i = 1 : size(barsOfColumns,1) 
    % 1 design individually bars of a column
    barNames = []; % to append in the end!!!
    for j = 1 : noStories %design of all bars of a column!
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesign(:,1,1) == barName);
        for k = 1 : length(cases)
            N_axial = DataDesign(barIndex, 2, k);
            My_h = DataDesign(barIndex, 3, k);
            Mz_b = DataDesign(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b); countf = countf +1;%, width);%, longRebarN, longRebarPh, givenLong)
            mAux1(k,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
        end
        [~, index] = max(mAux1(:,6));%best iteration for that bar
        bestIndi(j,:) = mAux1(index,:);%best individual bars of that column
    end
%
% the snippet
%     % 2 for each bar, design the column according to that one, each iteration is a page(3D)
    for k = 1 : noStories
        width = bestIndi(k,1); %let's call it parent bar
        longRebarN = bestIndi(k,3);
        longRebarPh = bestIndi(k,4);
        givenLong = columnComp(bestIndi(k,5));
        
        for j = 1 : noStories
            barName = barsOfColumns(i,j);
            barIndex = find(DataDesign(:,1,1) == barName);
            for p = 1 : length(cases)
                N_axial = DataDesign(barIndex, 2, p);
                My_h = DataDesign(barIndex, 3, p);
                Mz_b = DataDesign(barIndex, 4, p);
                [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, width, givenLong); countf = countf +1;%longRebarN, longRebarPh, 
                ratio = bestIndi(j,6) / reinfPercFin;
                mAux2(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, ratio, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
            end
            [~, index] = max(mAux2(:,6));%best iteration for that bar
            bestIndi2(j,:,k) = mAux2(index,:);%best individual bars of that column with a seed provided by the parent bar
        end
    end
%
    % 3 design from bottom to top
    %bottom bar
    bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration    %gWidth = floor(bigOrigWidth * .9 * 20)/20;
    gWidth = bigOrigWidth;
    barName = barsOfColumns(i,1);
    barIndex = find(DataDesign(:,1,1) == barName);
    for p = 1 : length(cases)
        N_axial = DataDesign(barIndex, 2, p);
        My_h = DataDesign(barIndex, 3, p);
        Mz_b = DataDesign(barIndex, 4, p);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth); countf = countf +1;%, longRebarN, longRebarPh, givenLong);
        ratio = bestIndi(1,6) / reinfPercFin;
        mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];%mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, ratio, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
    end
    [~, index] = max(mAux3(:,6));%best iteration for that bar
    bestIndi2(1,:,k+1) = mAux3(index,:);
            
    %other bars
    gWidth = floor(mAux3(index,1) * .9 * 20)/20;    %givenLong = mAux3(index, 5);
    for j = 2 : noStories
        barName = barsOfColumns(i,j);
        barIndex = find(DataDesign(:,1,1) == barName);
        for p = 1 : length(cases)
            N_axial = DataDesign(barIndex, 2, p);
            My_h = DataDesign(barIndex, 3, p);
            Mz_b = DataDesign(barIndex, 4, p);
            givenLong = columnComp(mAux3(index, 5),'EC8');
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong); countf = countf +1;%longRebarN, longRebarPh, 
            ratio = bestIndi(j,6) / reinfPercFin;
            mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, ratio, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
        end
        [~, index] = max(mAux3(:,6));%best iteration for that bar
        bestIndi2(j,:,k+1) = mAux3(index,:);%best individual bars of that column with a seed provided by the parent bar
        gWidth = floor(mAux3(index,1) * .9 * 20)/20;
        givenLong = mAux3(index, 6);
    end
    
    possDesigns = [];
    for z = 1 : (noStories + 1) % check if superior dimesions are the same or smaller, if not, kick it
        agr = issorted(bestIndi2(:,1,z),'descend');
        if agr == 1
            possDesigns = [possDesigns, z];
        end
    end
    validDesigns = bestIndi2(:,:,possDesigns);
        
    for z = 1 : size(validDesigns,3)
        eval(z) = mean(bestIndi2(:,7,z));
    end
    [~, bestDesign]= min(eval);
    
    columns = [columns; [barNames, validDesigns(:, [1:6,8:12],bestDesign)]];
    
    count = count + 1;
%     disp('');
%     disp(['Finished the ' num2str(count) ' column out of ' num2str(size(barsOfColumns,1))]);
%     disp(['Calls of beamDesign function: ' num2str(countf)]);
    toc
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress');% pause(1);
    if count == 2; break; end
end

for i = 1 : size(nodes,1)
    if nodes(i, 5) == 0 | nodes(i, 5) == noStories
        continue
    end
    
    %beams 

    
end

%toSeismo(columns, beams, nodes, element, stories, fck, fyk, cover)
close(loading);
disp('Finished');
toc