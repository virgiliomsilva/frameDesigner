%% DC3frameDesigner
function [] = DC3frameDesignerWSlab(buildingName, fck, fyk, cover, seismicVerticalLoadCase, seismicCases, nonSeismicCases, folder, flag, slabTopReinf)
%factor accounting for overstrength due to steel strain hardening and
%confinement of the concrete of the compression zone of the section
G_RD = 1.1;

tic; disp('Started DC3 slab');
%% STEP 1 - GETTING DATA
loading = waitbar(0,'Reading data','Name', 'DC3Slab: Step 1 of 9');

fnData = ['data\' buildingName '\dataset.csv'] ;
fnNodes = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesignMax, DataDesignMin, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);
clear buildingName fnData fnElement fnNodes

allCasesIdx = [1:length(cases)];
seismicVerticalLoadCaseIdx = find(ismember(cases, seismicVerticalLoadCase));
seismicCasesIdx = find(ismember(cases, seismicCases));
nonSeismicCasesIdx = find(ismember(cases, nonSeismicCases));

loading = waitbar(1,loading,'Reading data','Name', 'DC3Slab: Step 1 of 9'); pause(.5);
%% STEP 2 - DESIGN BEAMS
close(loading); loading = waitbar(0,'Initializing beams','Name', 'DC3Slab: Step 2 of 9'); pause(1);
beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesignMax(:,1,1) == beamDesiOrd(i));
    %longitudinal rebar
    mAux = [];
    for j = allCasesIdx
        M_Ed = DataDesignMax(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC3beamDesign(fck, fyk , cover, M_Ed, 0);
        mAux = [mAux;[sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition]];
        
        M_Ed = DataDesignMin(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC3beamDesign(fck, fyk , cover, M_Ed, 0);
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
    
    [~, ~, ~, ~, ~, ~, ~, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition] = DC3beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
    %midStirrups
    [shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid] = DC3beamDesignMidShear(fck, fyk , cover, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
    
    beams(end+1,:) = [DataDesignMax(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, mAux(conIndex, 7), shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition, Fz_Ed];
    beamsMid(end+1,:) = [DataDesignMax(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid, Fz_Ed];
    
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'DC1: Step 2 of 4');
end
save([folder '\DC3beamsIt1.mat'],'beams');
save([folder '\DC3beamsIt1mid.mat'],'beamsMid'); clear beamsMid
%% STEP 3 - DESIGN COLUMNS
close(loading); loading = waitbar(0,'Initializing columns','Name', 'DC3Slab: Step 3 of 9'); pause(1);
noStories = max(stories(:,1)); columns = []; columnsMid = []; count = 0;
for i = 1 : size(barsOfColumns,1)
    % #1 design individually bars of a column
    barNames = []; % to append in the end!!!
    for j = 1 : noStories %design of all bars of a column!
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesignMax(:,1,1) == barName);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        
        mAux1 = [];
        for k = nonSeismicCasesIdx
            N_axial = DataDesignMax(barIndex, 2, k);
            My_h = DataDesignMax(barIndex, 5, k);
            Mz_b = DataDesignMax(barIndex, 6, k);
            V_Ed = DataDesignMax(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, k);
            My_h = DataDesignMin(barIndex, 5, k);
            Mz_b = DataDesignMin(barIndex, 6, k);
            V_Ed = DataDesignMin(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1 = [mAux1;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %choose best based on biggest dimension
        mAux1 = sortrows(mAux1,[5 1],{'descend' 'ascend'});
        [~, index] = max(mAux1(:,1));
        bestIndi(j,:) = mAux1(index,:);
    end
    clear areaRebar barIndex barName index M_Rd mAux1 minWidth My_h Mz_b N_axial ...
        noRebar phiRebar reinfPercFin sCondition sec_b sec_h shearReinfArea ...
        shearReinfLoops shearReinfPhi shearReinfSpac V_Ed V_Rd j k
    
    % #2 design from bottom to top - based on the biggest width of the best individuals
    %   bottom bar
    bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration
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
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        
        N_axial = DataDesignMin(barIndex, 2, p);
        My_h = DataDesignMin(barIndex, 5, p);
        Mz_b = DataDesignMin(barIndex, 6, p);
        V_Ed = DataDesignMin(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
        mAux3 = [mAux3;[sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
    end
    mAux3 = sortrows(mAux3,[6 1],{'descend' 'ascend'});
    [~, index] = max(mAux3(:,6));%best iteration for that bar
    auxColumns = mAux3(index,:);
    clear areaRebar barIndex barName bigOrigWidth gWidth M_Rd minWidth ...
        areaMinIt My_h Mz_b N_axial noRebar p phiRebar reinfPercFin sCondition sec_b sec_h ...
        givenLong shearReinfArea shearReinfLoops shearReinfPhi shearReinfSpac V_Ed V_Rd ...
        barName index mAux3 minWidth
    %   other bars
    
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
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
            
            N_axial = DataDesignMin(barIndex, 2, p);
            My_h = DataDesignMin(barIndex, 5, p);
            Mz_b = DataDesignMin(barIndex, 6, p);
            V_Ed = DataDesignMin(barIndex, 4, p);
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, bigOrigWidth);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3 = [mAux3; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed]];
        end
        %best iteration for that bar
        mAux3 = sortrows(mAux3,[6 1],{'descend' 'ascend'});
        [~, index] = max(mAux3(:,6));
        auxColumns(j,:) = mAux3(index,:);
    end
    clear areaRebar barIndex barName bestIndi bigOrigWidth gWidth j M_Rd mAux3 minWidth ...
        areaMinIt My_h Mz_b N_axial noRebar p phiRebar reinfPercFin sCondition sec_b sec_h ...
        givenLong shearReinfArea shearReinfLoops shearReinfPhi shearReinfSpac V_Ed V_Rd index
    
    % check if dimensions agree (no bigger dimension on top of the column)
    %   if it fails will try another method to design it
    sorted = issorted(auxColumns(:,1),'descend');
    if ~sorted
        disp('First design method failed')
        %         %         it will take a lot of time but will fix possibly it
        %         %         this method will design the column so many times as individual
        %         %         bars of that column and then, per page of the matrix will put a
        %         %         column designed based on each bar
        %         for k = 1 : noStories
        %             width = bestIndi(k,1); % parent bar
        %             givenLong = columnComp(bestIndi(k,5), 'EC8');
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
        %                     [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, Gwidth);%, givenLong);
        %                     ratio = bestIndi(j,6) / reinfPercFin;
        %                     mAux2(p,:) = [ratio, sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
        %                 end
        %                 %best iteration for that bar
        %                 [~, index] = max(mAux2(:,7));
        %                 %best individual bars of that column with a seed provided by the parent bar
        %                 bestIndiMethod(j,:,k) = mAux2(index,:);
        %             end
        %         end
        %         clear k width givenLong givenLong barName barIndex minWidth Gwidth N_axial ...
        %             My_h Mz_b V_Ed sec_h sec_b noRebar phiRebar areaRebar reinfPercFin M_Rd ...
        %             shearReinfPhi shearReinfSpac shearReinfLoops shearReinfArea V_Rd sCondition ...
        %             ratio mAux2 index j
        %
        %         % check if superior dimesions are the same or smaller, if OK, add it to possDesigns
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
        %         clear possDesigns aCheck z validDesigns eval
        for t = noStories : -1 : 2
            if auxColumns(t-1,2) < auxColumns(t,2)
                gWidth = auxColumns(t,2);
                areaMinIt = auxColumns(t-1, 5);
                givenLong = columnLongMin(areaMinIt);
                [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, gWidth, givenLong);
                auxColumns(t-1,[1:7]) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd];
            end
            clear b h areaRebar barID N_AxialN M_RdN mAuxN p
        end
    end
    
    for j = 1 : size(auxColumns,1)
        sec_h = auxColumns(j, 1);
        noRebar = auxColumns(j, 3);
        phiRebar = auxColumns(j, 4);
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid] = DC3columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar);
        auxColumnsMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid];
    end
    
    columns = [columns; [barNames, auxColumns]];
    columnsMid = [columnsMid; [barNames, auxColumns(:,[1:7]), auxColumnsMid]];
    
    clear auxColumns auxColumnsMid barNames bestIndi noRebar phiRebar ...
        sec_h shearReinfAreaMid shearReinfLoopsMid shearReinfPhiMid ...
        shearReinfSpacMid sorted j
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'DC3Slab: Step 3 of 9');
end

save([folder '\DC3columnsIt1.mat'],'columns');
save([folder '\DC3columnsIt1mid.mat'],'columnsMid');
% clear columnsMid count i
%%
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC3slab\DC3beamsIt1.mat')
% load('D:\Google Drive\Disco UPorto\MIECIVIL\5 Ano\2 Semestre\MatlabFiles\frameDesigner\output\regular_DC3slab\DC3columnsIt1.mat')

%%
%% STEP 4 - UPDATE COLUMNS RESISTING BENDING MOMENTS
% updating the MRD values
close(loading); loading = waitbar(0,'Initializing columns resisting bending moments update','Name', 'DC3Slab: Step 4 of 9'); pause(1);
for i = 1 : size(columns,1)
    b = columns(i, 3);
    h = columns(i, 2);
    areaRebar = columns(i, 6);
    
    mAux = [];
    for j = seismicCasesIdx
        N_Axial = DataDesignMax(DataDesignMax(:,1,1) == columns(i,1) , 2, j);
        M_Rd = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
        mAux(end+1) = M_Rd;
    end
    
    columns(i, 8) = min(mAux);
    loading = waitbar(i / size(columns,1),loading,'Columns resisting bending moments update progress','Name', 'DC3Slab: Step 4 of 9');
end
clear areaRebar b h i j M_Rd mAux N_Axial
%% STEP 5 - DO THE RESISTING BENDING MOMENTS [1.3] COMPARISON
close(loading); loading = waitbar(0,'Initializing bending comparisons','Name', 'DC3Slab: Step 5 of 9'); pause(1);
minX = min(nodes(:,2)); maxX = max(nodes(:,2));
minY = min(nodes(:,3)); maxY = max(nodes(:,3));
borders = [minX,maxX,minY,maxY];

increNeed = [];
for i = 1 : size(nodes,1)
    %skip base points and on the top floor the verification is not needed
    if nodes(i, 5) == 0 %|| nodes(i, 5) == noStories
        continue
    end
    
    %get bars on that point
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 1); barsX = element(row,1);
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 2); barsY = element(row,1);
    [row, ~] = find(element(:,[2 3]) == nodes(i,1) & element(:,4) == 3); barsZ = element(row,1);
    
    %get design bending moments on those bars
    bendRdX = 0;
    for j = 1 : length(barsX)
        selfLength = DataDesignMax(DataDesignMax(:,1,1) == barsX(j) , 7, 1);
        plusLength = 0;
        for k = 1 : length(barsY)
            auxLength = DataDesignMax(DataDesignMax(:,1,1) == barsY(k) , 7, 1);
            plusLengthAux = min(selfLength/4, auxLength/2);
            plusLength = plusLength + plusLengthAux;
        end
        extraArea = plusLength * slabTopReinf;
        h = beams(beams(:,1) == barsX(j), 2);
        b = beams(beams(:,1) == barsX(j), 3);
        num = beams(beams(:,1) == barsX(j), 4);
        phi = beams(beams(:,1) == barsX(j), 5);
        A1 = pi * (phi/2000)^2 * num; A2 = A1;
        A1 = A1 + extraArea;
        
        MRD_beam = MrdBeam (fck, fyk, cover, b, h, A1, A2);
        bendRdX = bendRdX + MRD_beam;
    end
    
    bendRdY = 0;
    for j = 1 : length(barsY)
        selfLength = DataDesignMax(DataDesignMax(:,1,1) == barsY(j) , 7, 1);
        plusLength = 0;
        for k = 1 : length(barsX)
            auxLength = DataDesignMax(DataDesignMax(:,1,1) == barsX(k) , 7, 1);
            plusLengthAux = min(selfLength/4, auxLength/2);
            plusLength = plusLength + plusLengthAux;
        end
        extraArea = plusLength * slabTopReinf;
        h = beams(beams(:,1) == barsY(j), 2);
        b = beams(beams(:,1) == barsY(j), 3);
        num = beams(beams(:,1) == barsY(j), 4);
        phi = beams(beams(:,1) == barsY(j), 5);
        A1 = pi * (phi/2000)^2 * num; A2 = A1;
        A1 = A1 + extraArea;
        
        MRD_beam = MrdBeam (fck, fyk, cover, b, h, A1, A2);
        bendRdY = bendRdY + MRD_beam;
    end
    
    bendRdZ = 0;
    for j = 1 : length(barsZ)
        [columnRow, ~] = find(columns(:,1) == barsZ(j));
        bendRdZ = bendRdZ + columns(columnRow, 8);
    end
    
    %do the sums and comparisons
    ratioXZ = bendRdZ / bendRdX;
    if ratioXZ < 1.3
        auxBendMomX = bendRdX * 1.3;
    else
        auxBendMomX = 0;
    end
    
    ratioYZ = bendRdZ / bendRdY;
    if ratioYZ < 1.3
        auxBendMomY = bendRdY * 1.3;
    else
        auxBendMomY = 0;
    end
    increNeed(end+1,:) = [nodes(i,1), max(auxBendMomX, auxBendMomY)];
    
    clear auxBendMomX auxBendMomY barsX barsY barsZ beamRow bendRdX bendRdY ...
        bendRdZ columnRow j ratioXZ ratioYZ row
    
    loading = waitbar(i / (size(nodes,1)),loading,'Comparisons progress','Name', 'DC3Slab: Step 5 of 9');
end
save([folder '\increNeed.mat'],'increNeed');
save([folder '\columnsDel.mat'],'columns');
%% STEP 6 - REINFORCE COLUMNS FOR THE 1.3 RULE
close(loading); loading = waitbar(0,'Initializing column re-reinforcement','Name', 'DC3Slab: Step 6 of 9'); pause(1);
newColumns = []; columns13 = [];
for i = 1 : size(barsOfColumns,1)
    % build table
    aux1 = barsOfColumns(i,:);
    for j = 1 : length(aux1)
        auxElement(j,:) = element(element(:,1) == aux1(j),[1:3]);
    end
    
    %--NODE--%--BAR 1--%--MED 1--%--BAR 2--%--MED2--%--SUM MED--%--MED NEEDED--%--Diff--%
    table(:,1) = intersect([auxElement(:,2); auxElement(:,3)], increNeed(:,1));
    for j = 1 : size(table,1)
        [row, ~] = find(table(j,1) == auxElement(:, [2 3]));
        table(j, [2 4]) = auxElement(row ,1);
        table(j, 3) = columns(columns(:,1) == table(j, 2), 8);
        table(j, 5) = columns(columns(:,1) == table(j, 4), 8);
        table(j, 6) = table(j, 3) + table(j, 5);
        table(j, 7) = increNeed(table(j,1) == increNeed(:,1) , 2);
        table(j, 8) = table(j, 7) - table(j, 6);
    end
    
    clear auxElement
    
    % iterations
    while ~all(table(:,8) < 0)%while there is a "needed" bigger than what it is
        [~, index] = max(table(:, 8));
        toGive(:,1) = [table(index, 2); table(index,4)];
        toGive(:,2) = [table(index, 7)*(table(index, 3)/table(index, 6)); table(index, 7)*(table(index, 5)/table(index, 6))];
        
        %         improve = table(index, 7) / 2;
        %         if max(table(index, 3), table(index, 5)) > improve * 1.2 % if one if more than 60% of total needed just improve the weak one
        %             [~, col] = min([table(index, 3), table(index, 5)]);
        %             bar = table(index, col * 2);
        %             mImprove = table(index, 7) - max(table(index, 3), table(index, 5));
        %             toGive = [bar mImprove];
        %         else
        %             toGive(:,1) = [table(index, 2); table(index,4)];
        %             toGive(:,2) = [improve; improve];
        %         end
        
        % para a mesma taxa de armadura e sendo quadrados os pilares, a         % situação sísmica condicionante para o somatório é quando o        % esforço axial é máximo (das combinações sísmicas) pois reduz o
        % momendo resistente        % se eu dimensionar para o menor bending resistent        % necessário, só numa direção (hack: relação dos momentos tem de
        for j = 1 : size(toGive,1)
            matAux = [];
            for k = seismicCasesIdx
                N_axial = DataDesignMax(DataDesignMax(:,1,1) == toGive(j,1), 2, k);
                My_h = toGive(j, 2); Mz_b = 0;
                [columnsRow,~] = find(columns(:,1) == toGive(j,1));
                gWidth = columns(columnsRow, 2);
                areaOld = columns(columnsRow, 6);
                givLong = columnLongMin(areaOld);
                
                [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givLong);
                
                V_Ed = columns(columnsRow, 15);
                matAux = [matAux;toGive(j, 1), sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
            end
            matAux = sortrows(matAux,[6 2],{'descend' 'ascend'});
            [~, indexu] = max(matAux(:, 2));
            finalColumn = matAux(indexu,:);
            newColumns = [newColumns; finalColumn];
            
            %update my table!
            [updateRow, updateCol] = find(table(:,[2,4]) == toGive(j,1));
            updateCol = updateCol * 2;
            
            for k = 1 : length(updateCol)
                table(updateRow(k), [updateCol(k), updateCol(k) + 1]) = [toGive(j, 1), max(matAux(:, 8))];
            end
            
            table(:, 6) = table(:, 3) + table(:, 5);
            table(:, 8) = table(:, 7) - table(:, 6);
        end
        
        clear toGive areaOld areaRebar finalColumn givLong gWidth improve ...
            index indexu k longReinforce M_Rd My_h Mz_b N_axial noRebar ...
            phiRebar reinfPercFin row sCondition sec_b sec_h col mImprove ...
            shearReinfArea shearReinfLoops shearReinfPhi shearReinfSpac ...
            updateCol updateRow V_Ed V_Rd columnsRow
        
    end
    
    clear matAux
    
    %add the ones without changes and clean extra bars
    if ~isempty(newColumns)
        nonReReinf = setxor(aux1,newColumns(:,1));
        for j = 1 : length(nonReReinf)
            mAux = columns(columns(:,1) == nonReReinf(j),:);
            newColumns = [newColumns; mAux];
        end
        
        for j = 1 : noStories
            lastIndex = find(aux1(j) == newColumns(:,1), 1, 'last');
            primeColumns(j,:) = newColumns(lastIndex,:);
        end
    else
        primeColumns = newColumns;
    end
    clear auxElement aux1 table toGive newColumns lastIndex mAux nonReReinf
    
    %give concrete dimensions
    for z = 1 : (noStories) % check if superior dimesions are the same or smaller, if is it, add it to possDesigns
        aCheck = issorted(primeColumns(:,2),'descend');
        if ~aCheck
            for t = noStories : -1 : 2
                if primeColumns(t-1,2) < primeColumns(t,2)
                    gWidth = primeColumns(t,2);
                    areaMinIt = primeColumns(t-1, 6);
                    givenLong = columnLongMin(areaMinIt);
                    [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, gWidth, givenLong);
                    primeColumns(t-1,[2:8]) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd];
                end
            end
        end
    end
    
    clear aCheck z
    
    %   used to compatibilize reinforcement patterns allong height
    if flag == "yes"
        %for the given areas give the reinforcements
        %all the reinforcement possibilities
        auxVals = primeColumns(:, [1,2,6]);
        auxRefs = [];
        for j = 1 : noStories
            switch j
                case 1
                    auxRefs(:, [1:3], j) = auxVals;
                    width = auxVals(1,2);
                    reinfArea = auxVals(1,3);
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);
                    auxRefs(1, [4:6], j) = [phiRebar, noRebar, areaRebar];
                    givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    for k = 2 : noStories
                        width = auxVals(k,2);
                        reinfArea = auxVals(k,3);
                        [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea, givenLong);
                        auxRefs(k, [4:6], j) = [phiRebar, noRebar, areaRebar];
                        givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    end
                    
                case noStories
                    auxRefs(:, [1:3], j) = auxVals;
                    width = auxVals(noStories,2);
                    reinfArea = auxVals(noStories,3);
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);
                    auxRefs(noStories, [4:6], j) = [phiRebar, noRebar, areaRebar];
                    givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    for k = (noStories - 1) : -1 : 1
                        width = auxVals(k,2);
                        reinfArea = auxVals(k,3);
                        [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea, givenLong);
                        auxRefs(k, [4:6], j) = [phiRebar, noRebar, areaRebar];
                        givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    end
                    
                otherwise
                    auxRefs(:, [1:3], j) = auxVals;
                    width = auxVals(j,2);
                    reinfArea = auxVals(j,3);
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);
                    auxRefs(j, [4:6], j) = [phiRebar, noRebar, areaRebar];
                    givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    for k = (j+1) : noStories
                        width = auxVals(k,2);
                        reinfArea = auxVals(k,3);
                        [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea, givenLong);
                        auxRefs(k, [4:6], j) = [phiRebar, noRebar, areaRebar];
                        givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    end
                    
                    width = auxVals(j,2);
                    reinfArea = auxVals(j,3);
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);
                    auxRefs(j, [4:6], j) = [phiRebar, noRebar, areaRebar];
                    givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    for k = (j-1) : -1 : 1
                        width = auxVals(k,2);
                        reinfArea = auxVals(k,3);
                        [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea, givenLong);
                        auxRefs(k, [4:6], j) = [phiRebar, noRebar, areaRebar];
                        givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                    end
            end
            
            clear areaRebar auxVals j k noRebar phiRebar reinfArea width
        end
        
        %delete the impossoble ones
        [~, delIdx] = find(auxRefs(:,6,:) == 0);
        if length(unique(delIdx)) == noStories
            for k = 1 : size(auxRefs,3)
                for j = 1 : size(auxRefs,1)
                    if auxRefs(j,4,k) == 0
                        width = auxRefs(j,2,k);
                        reinfArea = auxRefs(j,3,k);
                        [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);
                        auxRefs(j, [4:6], k) = [phiRebar, noRebar, areaRebar];
                    end
                end
            end
            delIdx = [];
        else
            auxRefs(:, :, delIdx) = [];
        end
        
        clear delIdx
        
        %calculate the overshooting and the best in average
        auxRefs(:,7,:) = auxRefs(:,6,:) ./ auxRefs(:,3,:);
        means = mean(auxRefs(:,7,:),1);
        [~, bestIdx] = min(means);
        auxRefs = auxRefs(:,:,bestIdx);
        
        %update columns
        for j = 1 : noStories
            primeColumns(j, [4:6]) = auxRefs(j, [4:6]);
            areaRebar = primeColumns(j, 6);
            sec_b = primeColumns(j,2); sec_h = sec_b;
            primeColumns(j, 7) = areaRebar * fyk / (sec_b  * sec_h * fck);
            N_axial = max(DataDesignMax(DataDesignMax(:,1,1) == primeColumns(j,1), 2, seismicCasesIdx));
            primeColumns(j, 8) = MrdColumn(fck, fyk, sec_b, sec_h, areaRebar, N_axial);
            
            longReinfN = primeColumns(j, 5);
            longReinfPh = primeColumns(j, 4);
            Ved = columns(columns(:,1) == primeColumns(j, 1), 15);
            [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, Vrd, sCondition] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, sec_b, [longReinfPh,longReinfN,areaRebar], [], sec_b, longReinfN, longReinfPh);
            primeColumns(j, [9:15]) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, Vrd, sCondition, Ved];
        end
        
        clear auxRefs areaRebar bestIdx givenLong longReinfN longReinfPh means ...
            N_axial sCondition sec_b sec_h shearReinfArea shearReinfLoops ...
            shearReinfPhi shearReinfSpac Ved Vrd
    end
    
    columns13 = [columns13; primeColumns];
    
    loading = waitbar(i / (size(barsOfColumns,1)),loading,'Column re-reinforcement progress','Name', 'DC3Slab: Step 6 of 9');
    clear primeColumns
    newColumns = [];
end
save([folder '\DC3columnsIt2rule13.mat'],'columns13');
clear noStories flag columns increNeed newColumns
%% STEP 7 - SHEAR CAPACITY DESIGN BEAMS
close(loading); loading = waitbar(0,'Initializing Beam seismic equilibrium','Name', 'DC3Slab: Step 7 of 9'); pause(1);
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
                [barRow, ~] = find(columns13(:,1) == barsZ(j));
                b = columns13(barRow, 3);
                h = columns13(barRow, 2);
                areaRebar = columns13(barRow, 6);
                N_Axial = DataDesignMax(DataDesignMax(:,1) == barsZ(j) , 2, p);
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
    Vshear = max(max(Vseismo));
    
    % shear reinforcement
    barIndex = find(beams(:,1) == beamDesiOrd(i));
    given_b = beams(barIndex, 3);
    given_h = beams(barIndex, 2);
    longReinfN =  beams(barIndex, 4);
    longReinfPh =  beams(barIndex, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition] = DC3beamDesign(fck, fyk , cover, 10, Vshear, given_b, given_h, longReinfN, longReinfPh);
    seismicBeams = [seismicBeams; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition, Vshear]];
    
    %mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = DC3beamDesignMidShear(fck, fyk , cover, Vshear, given_b, given_h, longReinfN, longReinfPh);
    seismicBeamsMidShear = [seismicBeamsMidShear; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, 0, Vshear]];
    
    clear barIndex cLength given_b given_h longReinfN longReinfPh Mid sCondition ...
        shearReinfLoops shearReinfPhi shearReinfSpac V_Rd Vseismic Vshear
    
    loading = waitbar(i / (length(beamDesiOrd)),loading,'Beam seismic equilibrium','Name', 'DC3Slab: Step 7 of 9');
end
save([folder '\DC3beamsIt2SeismicEq.mat'],'seismicBeams');
save([folder '\DC3beamsIt2SeismicEqMid.mat'],'seismicBeamsMidShear');
clear seismicBeamsMidShear
%% STEP 8 - SHEAR CAPACITY DESIGN COLUMNS
close(loading); loading = waitbar(0,'Initializing column seismic equilibrium','Name', 'DC3Slab: Step 8 of 9'); pause(1);
pilars = setxor(beamDesiOrd, element(:,1));
seismicColumns = []; seismicColumnsMidShear = [];
for i = 1 : length(pilars)
    %bar info
    barRow = find(columns13(:,1) == pilars(i));
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
        b = columns13(barRow, 3);
        h = columns13(barRow, 2);
        areaRebar = columns13(barRow, 6);
        N_Axial = DataDesignMax(DataDesignMax(:,1) == pilars(i) , 2, p);
        Mrc = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
        
        %obtain max sum of bending moment on top each node of the column
        for k = 1 : 2
            [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
            bendRdZaux = 0;
            for j = 1 : length(barsZ)
                b = columns13(columns13(:,1) == barsZ(j), 3);
                h = columns13(columns13(:,1) == barsZ(j), 2);
                areaRebar = columns13(columns13(:,1) == barsZ(j), 6);
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
    Vshear = max(max(Vseismo));
    
    given_h = columns13(barRow, 2);
    longReinfN = columns13(barRow, 4);
    longReinfPh = columns13(barRow, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, given_h, [longReinfPh,longReinfN,areaRebar], Vshear, given_h, longReinfN, longReinfPh);
    seismicColumns = [seismicColumns; [columns13(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, Vshear]];
    
    %seismic beams mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC3columnDesignMidShear(fck, fyk , cover, given_h, longReinfN, longReinfPh);
    seismicColumnsMidShear = [seismicColumnsMidShear; [columns13(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, 0, 0, 0]];
    
    clear areaRebar b barRow bendMomsaux bendMomsBeam bendRdX bendRdY ...
        bendRdZ cLength h j k longReinfN longReinfPh M_RdS mAuxS maxIndex Mid ...
        MRc N_AxialS p shearReinfArea shearReinfLoops shearReinfPhi shearReinfSpac ...
        V_Rd Vshear
    
    loading = waitbar(i / (length(pilars)),loading,'Column seismic equilibrium progress','Name', 'DC3Slab: Step 8 of 9');
end
save([folder '\DC3columnsIt3SeismicEq.mat'],'seismicColumns');
save([folder '\DC3columnsIt3SeismicEqMid.mat'],'seismicColumnsMidShear');
%% STEP 9 - EXPORT RESULTS TO CSV FOR SEISMOSTRUCT
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'DC3Slab: Step 9 of 9'); pause(1);
toSeismo(seismicColumns(:,[1:5, 9:11]), seismicBeams(:,[1:5, 8:10]), nodes, element, stories, fck, fyk, cover, folder);
loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'DC3Slab: Step 9 of 9');
%% FINISH
time = toc; save([folder '\time.mat'],'time');
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);