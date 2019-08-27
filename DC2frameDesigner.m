%% DC2frameDesigner
function [] = DC2frameDesigner(buildingName, fck, fyk, cover, seismicVerticalLoadCase, seismicCases, nonSeismicCases, folder)
%factor accounting for overstrength due to steel strain hardening and 
%confinement of the concrete of the compression zone of the section
G_RD = 1.1;

tic; disp('Started DC2');
loading = waitbar(0,'Reading data','Name', 'DC2: Step 1 of 6');

fnData = ['data\' buildingName '\dataset.csv'] ;
fnNodes = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesignMax, DataDesignMin, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);

allCasesIdx = [1:length(cases)];
seismicVerticalLoadCaseIdx = find(ismember(cases, seismicVerticalLoadCase));
seismicCasesIdx = find(ismember(cases, seismicCases));
nonSeismicCasesIdx = find(ismember(cases, nonSeismicCases));

loading = waitbar(1,loading,'Reading data','Name', 'DC2: Step 1 of 6'); pause(.5);
%%
close(loading); loading = waitbar(0,'Initializing beams','Name', 'DC2: Step 2 of 6'); pause(1);
beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesignMax(:,1,1) == beamDesiOrd(i));
    
    %longitudinal rebar
    mAux = [];
    for j = allCasesIdx
        M_Ed = DataDesignMax(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC2beamDesign(fck, fyk , cover, M_Ed, 0);
        mAux = [mAux;[sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition]];
        
        M_Ed = DataDesignMin(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC2beamDesign(fck, fyk , cover, M_Ed, 0);
        mAux = [mAux;[sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition]];
    end
    mAux = sortrows(mAux, [2 1 5]);
    [M_rd, conIndex] = max(mAux(:,6)); %best M_rd
    
    %stirrups
    given_h = mAux(conIndex, 1);
    given_b = mAux(conIndex, 2);
    longRebarN = mAux(conIndex, 3);
    longRebarPh = mAux(conIndex, 4);
    
    Fz_Ed = max(max(abs(DataDesignMax(barIndex, 4, allCasesIdx)),abs(DataDesignMin(barIndex, 4, allCasesIdx))));
    
    [~, ~, ~, ~, ~, ~, ~, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition] = DC2beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
    %midStirrups
    [shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid] = DC2beamDesignMidShear(fck, fyk , cover, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);

    beams(end+1,:) = [DataDesignMax(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, mAux(conIndex, 7), shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition, Fz_Ed];
    beamsMid(end+1,:) = [DataDesignMax(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid, Fz_Ed];
    
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'DC1: Step 2 of 4');
end
save([folder '\DC2beamsIt1.mat'],'beams')
save([folder '\DC2beamsIt1mid.mat'],'beamsMid')
%%
close(loading); loading = waitbar(0,'Initializing columns','Name', 'DC2: Step 3 of 6'); pause(1);
noStories = max(stories(:,1)); count = 0; columns = []; columnsMid = [];
for i = 1 : size(barsOfColumns,1)
    % #1 design individually bars of a column
    barNames = []; % to append in the end
    
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
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesignIts(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, k);
            My_h = DataDesignMin(barIndex, 5, k);
            Mz_b = DataDesignMin(barIndex, 6, k);
            V_Ed = DataDesignMin(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesignIts(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %choose best based on biggest dimension
        mAux1 = sortrows(mAux1,[5 1],{'descend' 'ascend'});
        [~, index] = max(mAux1(:,1));
        bestIndi(j,:) = mAux1(index,:);
    end
    
    % #2 design from bottom to top - based on the biggest width of the best individuals
    %   bottom bar
    bigOrigWidth = max(bestIndi(:,1));
    barName = barsOfColumns(i,1);
    barIndex = find(DataDesignMax(:,1,1) == barName);
    areaMinIt = bestIndi(1, 5);
    givenLong = columnLongMin(areaMinIt);
    mAux3 = [];
    for p = allCasesIdx
        N_axial = DataDesignMax(barIndex, 2, p);
        My_h = DataDesignMax(barIndex, 5, p);
        Mz_b = DataDesignMax(barIndex, 6, p);
        V_Ed = DataDesignMax(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        
        N_axial = DataDesignMin(barIndex, 2, p);
        My_h = DataDesignMin(barIndex, 5, p);
        Mz_b = DataDesignMin(barIndex, 6, p);
        V_Ed = DataDesignMin(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
    end
    
    %best iteration for that bar
    mAux3 = sortrows(mAux3,[5 1],{'descend' 'ascend'});
    [~, index] = max(mAux3(:,6));
    auxColumns = mAux3(index,:);
    
    %   other bars
%     try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
%     gWidth = max(minWidth, min(floor(mAux3(index,1) * .9 * 20)/20, mAux3(index,1) - .05));
    
    for j = 2 : noStories
        bigOrigWidth = bestIndi(j,1);
        barName = barsOfColumns(i,j);
        barIndex = find(DataDesignMax(:,1,1) == barName);
        areaMinIt = bestIndi(j, 5);
        givenLong = columnLongMin(areaMinIt);
        mAux3 = [];
        for p = allCasesIdx
            N_axial = DataDesignMax(barIndex, 2, p);
            My_h = DataDesignMax(barIndex, 5, p);
            Mz_b = DataDesignMax(barIndex, 6, p);
            V_Ed = DataDesignMax(barIndex, 4, p);
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, bigOrigWidth);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, p);
            My_h = DataDesignMin(barIndex, 5, p);
            Mz_b = DataDesignMin(barIndex, 6, p);
            V_Ed = DataDesignMin(barIndex, 4, p);
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, bigOrigWidth);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %best iteration for that bar
        mAux3 = sortrows(mAux3,[5 1],{'descend' 'ascend'});
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
        %         %         this method will design the column so many times as individual
        %         %         bars of that column and then, per page of the matrix will put a
        %         %         column designed based on each bar
        %         for k = 1 : noStories
        %             % parent bar
        %             width = bestIndi(k,1);
        % %             givenLong = columnComp(bestIndi(k,5), 'EC8');
        %
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
        %                     [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, Gwidth);
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
            if auxColumns(t-1,2) < auxColumns(t,2)
                gWidth = auxColumns(t,2);
                areaMinIt = auxColumns(t-1, 5);  
                givenLong = columnLongMin(areaMinIt);
                [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC2columnDesign(fck, fyk , cover, 10, 10, 10, gWidth, givenLong);
                auxColumns(t-1,[1:7]) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd];
            end
        end
    end
    
    for j = 1 : size(auxColumns,1)
        sec_h = auxColumns(j, 1);
        noRebar = auxColumns(j, 3);
        phiRebar = auxColumns(j, 4);
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid] = DC2columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar);
        auxColumnsMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid];
    end
    
    columns = [columns; [barNames, auxColumns]];
    columnsMid = [columnsMid; [barNames, auxColumns(:,[1:7]), auxColumnsMid]];
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'DC2: Step 3 of 6');
end
save([folder '\DC2columnsIt1.mat'],'columns')
save([folder '\DC2columnsIt1mid.mat'],'columnsMid')
%%
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC2\DC2columnsIt1.mat')
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC2\DC2columnsIt1mid.mat')
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC2\DC2beamsIt1mid.mat')
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC2\DC2beamsIt1.mat')
%% seismic equilibrium
%BEAM
close(loading); loading = waitbar(0,'Initializing beam seismic equilibrium','Name', 'DC2: Step 4 of 6'); pause(1);
seismicBeams = []; seismicBeamsMidShear = [];
for i = 1 : length(beamDesiOrd)
    %bar info
    res = element(element(:,1) == beamDesiOrd(i), [2 3]);
    node1 = res(1); node2 = res(2);
    auxNode = [node1, node2];
    cLength = DataDesignMax(DataDesignMax(:,1,1) == beamDesiOrd(i) , 7, 1);
    direction = element(element(:,1) == beamDesiOrd(i), 4);
    
    for k = 1 : 2
        %obtain max sum of bending moment on top each node of the beams, each direction
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == direction); barsBeam = element(row,1);
        bendRdNodeAux = 0;
        for j = 1 : length(barsBeam)
            M_Rd = beams(beams(:,1) == barsBeam(j), 6);
            bendRdNodeAux = bendRdNodeAux + M_Rd;
        end
        bendRdBeam(k) = bendRdNodeAux;
    end
    
    % for each combination
    Vseismo = [];
    for p = seismicCasesIdx
        for k = 1 : 2
            % for each node get the bars that connect there
            [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
            % get the sumMrd for each node
            bendRdZAux = 0;
            for j = 1 : length(barsZ)
                [barRow, ~] = find(columns(:,1) == barsZ(j));
                b = columns(barRow, 3);
                h = columns(barRow, 2);
                areaRebar = columns(barRow, 6);
                N_Axial = DataDesignMax(barRow , 2, p);
                Mrc = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
                bendRdZAux = bendRdZAux + Mrc;
            end
            bendRdZ(k) = bendRdZAux;
            
            if bendRdBeam(k) < bendRdZ(k)
                Mid(k) = G_RD * bendRdBeam(k);
            else
                Mid(k) = G_RD * bendRdBeam(k) * (bendRdZ(k) / bendRdBeam(k));
            end
        end
        
        VofMom = sum(Mid)/cLength;
        V1 = DataDesignMax(DataDesignMax(:,1,1) == beamDesiOrd(i), 4, seismicVerticalLoadCaseIdx);
        V2 = DataDesignMin(DataDesignMin(:,1,1) == beamDesiOrd(i), 4, seismicVerticalLoadCaseIdx);
        Vseismo = [Vseismo; V1 + VofMom];
        Vseismo = [Vseismo; V2 - VofMom];
    end
        
    Vseismo = abs(Vseismo);
    Vshear = max(Vseismo);
    
    % shear reinforcement
    barIndex = find(beams(:,1) == beamDesiOrd(i));
    given_b = beams(barIndex, 3);
    given_h = beams(barIndex, 2);
    longReinfN =  beams(barIndex, 4);
    longReinfPh =  beams(barIndex, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition] = DC2beamDesign(fck, fyk , cover, 10, Vshear, given_b, given_h, longReinfN, longReinfPh);
    seismicBeams = [seismicBeams; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition, Vshear]];
    
    %seismic beams mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = DC2beamDesignMidShear(fck, fyk , cover, Vshear, given_b, given_h, longReinfN, longReinfPh);
    seismicBeamsMidShear = [seismicBeamsMidShear; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, 0, Vshear]];
        
    loading = waitbar(i / (length(beamDesiOrd)),loading,'Beam seismic equilibrium','Name', 'DC2: Step 4 of 6');
end
save([folder '\DC2beamsIt2SeismicEq.mat'],'seismicBeams')
save([folder '\DC2beamsIt2SeismicEqMid.mat'],'seismicBeamsMidShear')

%COLUMN
close(loading); loading = waitbar(0,'Initializing column seismic equilibrium','Name', 'DC2: Step 5 of 6'); pause(1);
pilars = setxor(beamDesiOrd, element(:,1));
seismicColumns = []; seismicColumnsMidShear = [];
for i = 1 : length(pilars)
    %bar info
    barRow = find(columns(:,1) == pilars(i));
    res = element(element(:,1) == pilars(i), [2 3]);
    node1 = res(1); node2 = res(2);
    auxNode = [node1, node2];
    cLength = DataDesignMax(DataDesignMax(:,1,1) == pilars(i) , 7, 1);
    
    for k = 1 : 2
        %obtain max sum of bending moment on top each node of the beams, each direction
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 1); barsX = element(row,1);
        bendRdXaux = 0;
        for j = 1 : length(barsX)
            M_Rd = beams(beams(:,1) == barsX(j), 6);
            bendRdXaux = bendRdXaux + M_Rd;
        end
        bendRdX(k) = bendRdXaux;
        
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 2); barsY = element(row,1);
        bendRdYaux = 0;
        for j = 1 : length(barsY)
            M_Rd = beams(beams(:,1) == barsY(j), 6);
            bendRdYaux = bendRdYaux + M_Rd;
        end
        bendRdY(k) = bendRdYaux;
    end
       
    % for each combination
    Vseismo = [];
    for p = seismicCasesIdx
        %Mrc of the given bar
        b = columns(barRow, 3);
        h = columns(barRow, 2);
        areaRebar = columns(barRow, 6);
        N_Axial = DataDesignMax(barRow , 2, p);
        Mrc = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
        
        %obtain max sum of bending moment on top each node of the column
        for k = 1 : 2
            [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
            bendRdZaux = 0;
            for j = 1 : length(barsZ)
                b = columns(columns(:,1) == barsZ(j), 3);
                h = columns(columns(:,1) == barsZ(j), 2);
                areaRebar = columns(columns(:,1) == barsZ(j), 6);
                N_Axial = DataDesignMax(DataDesignMax(:,1,1) == barsZ(j) , 2, p);
                M_Rd = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
                bendRdZaux = bendRdZaux + M_Rd;
            end
            bendRdZ(k) = bendRdZaux;
        end
        
        for k = 1 : 2
            if bendRdX(k) > bendRdZ(k)
                MidX(k) = G_RD * Mrc;
            else
                MidX(k) = G_RD * Mrc * (bendRdX(k) / bendRdZ(k));
            end
        
            if bendRdY(k) > bendRdZ(k)
                MidY(k) = G_RD * Mrc;
            else
                MidY(k) = G_RD * Mrc * (bendRdY(k) / bendRdZ(k));
            end
        end
    Vseismo = [Vseismo; sum(MidX)/cLength];
    Vseismo = [Vseismo; sum(MidY)/cLength];
    end

    Vseismo = abs(Vseismo);
    Vshear = max(Vseismo);
    
    given_h = columns(barRow, 2);
    longReinfN = columns(barRow, 4);
    longReinfPh = columns(barRow, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, 10, 10, 10, given_h, [longReinfPh,longReinfN,areaRebar], Vshear, given_h, longReinfN, longReinfPh);
    seismicColumns = [seismicColumns; [columns(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, Vshear]];
    
    %seismic beams mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC2columnDesignMidShear(fck, fyk , cover, given_h, longReinfN, longReinfPh);
    seismicColumnsMidShear = [seismicColumnsMidShear; [columns(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, 0, 0, 0]];
    
    loading = waitbar(i / (length(pilars)),loading,'Column seismic equilibrium progress','Name', 'DC2: Step 5 of 6');
end
save([folder '\DC2columnsIt2SeismicEq.mat'],'seismicColumns')
save([folder '\DC2columnsIt2SeismicEqMid.mat'],'seismicColumnsMidShear')
%%
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'DC2: Step 6 of 6'); pause(1);
toSeismo(seismicColumns(:,[1:5, 9:11]), seismicBeams(:,[1:5, 8:10]), nodes, element, stories, DataDesignMin, DataDesignMax, seismicVerticalLoadCaseIdx, fck, fyk, cover, folder);
loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'DC2: Step 6 of 6');

time = toc; save([folder '\time.mat'],'time');
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);