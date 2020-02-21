%% DC1frameDesigner
function [] = DC1frameDesigner(buildingName, fck, fyk, cover, seismicVerticalLoadCase, nonSeismicCases, folder)
tic; disp('Started DC1');
loading = waitbar(0,'Reading data','Name', 'DC1: Step 1 of 4');

fnData    = ['data\' buildingName '\dataset.csv'] ;
fnNodes   = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesignMax, DataDesignMin, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);
clear buildingName fnData fnElement fnNodes

allCasesIdx = [1:length(cases)];
seismicVerticalLoadCaseIdx = find(ismember(cases, seismicVerticalLoadCase));
nonSeismicCasesIdx = find(ismember(cases, nonSeismicCases));

loading = waitbar(1,loading,'Reading data','Name', 'DC1: Step 1 of 4'); pause(1);
%%
close(loading); loading = waitbar(0,'Initializing beams','Name', 'DC1: Step 2 of 4'); pause(1);
beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesignMax(:,1,1) == beamDesiOrd(i));
    
    %longitudinal rebar and stirrups on critical region
    maxM_Ed = max(abs(DataDesignMax(barIndex, 5, :)));
    Fz_Ed = max(max(abs(DataDesignMax(barIndex, 4, allCasesIdx)),abs(DataDesignMin(barIndex, 4, allCasesIdx))));
    
    [sec_h, sec_b, longReinfNo, longReinfPhi, ~, M_Rd, roMinCondition, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition] = ...
        DC1beamDesign(fck, fyk , cover, maxM_Ed, Fz_Ed);
    
    %midStirrups
    [shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid] = DC1beamDesignMidShear(fck, fyk , cover, Fz_Ed, sec_b, sec_h, longReinfNo, longReinfPhi);
    
    beams   (end+1,:) = [DataDesignMax(barIndex,1,1), sec_h, sec_b, longReinfNo, longReinfPhi, M_Rd, roMinCondition, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition, Fz_Ed];
    beamsMid(end+1,:) = [DataDesignMax(barIndex,1,1), sec_h, sec_b, longReinfNo, longReinfPhi, M_Rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid, Fz_Ed];
    
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'DC1: Step 2 of 4');
end
save([folder '\DC1beamsIt1.mat'],   'beams')
save([folder '\DC1beamsIt1mid.mat'],'beamsMid')

clear barIndex beams beamsMid Fz_Ed i longReinfArea longReinfNo longReinfPhi M_Rd maxM_Ed roMinCondition sCondition ...
    sec_h sec_b shearLegs shearLegsMid shearPhi shearPhiMid shearSpac shearSpacMid V_Rd_it V_RdMid
%%
close(loading); loading = waitbar(0,'Initializing columns','Name', 'DC1: Step 3 of 4'); pause(1);
noStories = max(stories(:,1)); count = 0; columns = []; columnsMid = [];
for i = 1 : size(barsOfColumns,1)
    % #1 design individually bars of a column
    barNames = []; % to append in the end!!!
    
    %design of all bars of a column!
    for j = 1 : noStories
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesignMax(:,1,1) == barName);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        
        mAux1 = [];
        for k = nonSeismicCasesIdx
            N_axial = DataDesignMax(barIndex, 2, k);
            My_h = DataDesignMax(barIndex, 5, k);
            Mz_b = DataDesignMax(barIndex, 6, k);
            V_Ed = DataDesignMax(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesignIts(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, k);
            My_h = DataDesignMin(barIndex, 5, k);
            Mz_b = DataDesignMin(barIndex, 6, k);
            V_Ed = DataDesignMin(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesignIts(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %choose best based on biggest dimension
        mAux1 = sortrows(mAux1,[5 1],{'descend' 'ascend'});
        [~, index] = max(mAux1(:,1));
        bestIndi(j,:) = mAux1(index,:);
    end
    
    % #2 design from bottom to top - based on the biggest width and reinf
    % area of the best individuals - check if the chosen columns is good
    % for all the other combinations
    %   bottom bar
    bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration
    barName = barsOfColumns(i,1);
    barIndex = find(DataDesignMax(:,1,1) == barName);
    areaMinIt = bestIndi(1, 5);
    givenLong = columnLongMin(areaMinIt, 'flag');
    mAux3 = [];
    for p = allCasesIdx
        N_axial = DataDesignMax(barIndex, 2, p);
        My_h = DataDesignMax(barIndex, 5, p);
        Mz_b = DataDesignMax(barIndex, 6, p);
        V_Ed = DataDesignMax(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        
        N_axial = DataDesignMin(barIndex, 2, p);
        My_h = DataDesignMin(barIndex, 5, p);
        Mz_b = DataDesignMin(barIndex, 6, p);
        V_Ed = DataDesignMin(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
    end
    %best iteration for that bar, best reinforcement pattern
    mAux3 = sortrows(mAux3,[5 1],{'descend' 'ascend'});
    [~, index] = max(mAux3(:,6));
    auxColumns = mAux3(index,:);
    
    %   other bars
    for j = 2 : noStories
        bigOrigWidth = bestIndi(j,1);
        barName = barsOfColumns(i,j);
        barIndex = find(DataDesignMax(:,1,1) == barName);
        areaMinIt = bestIndi(j, 5);
        givenLong = columnLongMin(areaMinIt, 'flag');
        mAux3 = [];
        for p = allCasesIdx
            N_axial = DataDesignMax(barIndex, 2, p);
            My_h = DataDesignMax(barIndex, 5, p);
            Mz_b = DataDesignMax(barIndex, 6, p);
            V_Ed = DataDesignMax(barIndex, 4, p);
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, bigOrigWidth);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, p);
            My_h = DataDesignMin(barIndex, 5, p);
            Mz_b = DataDesignMin(barIndex, 6, p);
            V_Ed = DataDesignMin(barIndex, 4, p);
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, bigOrigWidth);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %best iteration for that bar
        mAux3 = sortrows(mAux3,[6 1],{'descend' 'ascend'});
        [~, index] = max(mAux3(:,6));
        auxColumns(j,:) = mAux3(index,:);
        
        %         try
        %             try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        %             gWidth = max(minWidth, min(floor(mAux3(index,1) * .9 * 20)/20, mAux3(index,1) - .05));
        %             %             givenLong = mAux3(index, 6);
        %         end
    end
    
    % check if dimensions agree (no bigger dimension on top of the column)
    % if it fails will try another method to design it
    sorted = issorted(auxColumns(:,1),'descend');
    if ~sorted
        disp('First design method failed')
        %         this method will design the column so many times as individual
        %         bars of that column and then, per page of the matrix will put a
        %         column designed based on each bar
        %         for k = 1 : noStories
        %             % parent bar
        %             width = bestIndi(k,1);
        %             %             givenLong = columnComp(bestIndi(k,5), 'EC8');
        %             for j = 1 : noStories
        %                 barName = barsOfColumns(i,j);
        %                 barIndex = find(DataDesign(:,1,1) == barName);
        %                 try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        %                 Gwidth = max(width, minWidth);
        %                 for p = 1 : length(cases)
        %                     N_axial = DataDesign(barIndex, 2, p);
        %                     My_h = DataDesign(barIndex, 5, p);
        %                     Mz_b = DataDesign(barIndex, 6, p);
        %                     V_Ed = DataDesign(barIndex, 4, p);
        %                     [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, Gwidth);%, givenLong);
        %                     ratio = bestIndi(j,6) / reinfPercFin;
        %                     mAux2(p,:) = [ratio, sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
        %                 end
        %                 %best iteration for that bar
        %                 [~, index] = max(mAux2(:,7));
        %                 %best individual bars of that column with a seed provided by the parent bar
        %                 bestIndiMethod(j,:,k) = mAux2(index,:);
        %             end
        %         end
        %
        %         % check if superior dimesions are the same or smaller, if is it, add it to possDesigns
        %         possDesigns = [];
        %         for z = 1 : (noStories)
        %             aCheck = issorted(bestIndiMethod(:,2,z),'descend');
        %             if aCheck == 1
        %                 possDesigns = [possDesigns, z];
        %             end
        %         end
        %
        %         if isempty(possDesigns)
        %             disp(['Cannot design this bars/column: ' num2str(barNames')]);
        %             [row, col] = size(auxColumns);
        %             auxColumns = zeros(row, col);
        %         else
        %             validDesigns = bestIndiMethod(:,:,possDesigns);
        %             for z = 1 : size(validDesigns,3)
        %                 eval(z) = mean(bestIndiMethod(:,1,z));
        %             end
        %             [~, bestDesign]= min(eval);
        %             auxColumns = bestIndiMethod(:,[2:end],bestDesign);
        %         end
        
        for t = noStories : -1 : 2
            % give the needed dimension and check if minimum einforcement is provided
            if auxColumns(t-1,2) < auxColumns(t,2)
                gWidth = auxColumns(t,2);
                areaMinIt = auxColumns(t-1, 5);
                givenLong = columnLongMin(areaMinIt);
                [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC1columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
                auxColumns(t-1,[1:7]) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd];
            end
        end
    end
    
    for j = 1 : size(auxColumns,1)
        sec_h = auxColumns(j, 1);
        noRebar = auxColumns(j, 3);
        phiRebar = auxColumns(j, 4);
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid] = DC1columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar);
        auxColumnsMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid];
    end
    
    columns = [columns; [barNames, auxColumns]];
    columnsMid = [columnsMid; [barNames, auxColumns(:,[1:7]), auxColumnsMid]];
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'DC1: Step 3 of 4');
end
save([folder '\DC1columnsIt1.mat'],'columns')
save([folder '\DC1columnsIt1mid.mat'],'columnsMid')
%%
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'DC1: Step 4 of 4'); pause(1);
toSeismo(columns(:,[1:5, 9:11]), beams(:,[1:5, 8:10]), nodes, element, stories, DataDesignMin, DataDesignMax, seismicVerticalLoadCaseIdx, fck, fyk, cover, folder);
loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'DC1: Step 4 of 4');
%%
time = toc; save([folder '\time.mat'],'time');
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);