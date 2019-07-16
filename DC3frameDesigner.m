%% DC3frameDesigner
function [] = DC3frameDesigner(buildingName, fck, fyk, cover, seismicCases, folder)
%factor accounting for overstrength due to steel strain hardening and 
%confinement of the concrete of the compression zone of the section
G_RD = 1.1; 
tic; disp('Started');
loading = waitbar(0,'Reading data','Name', 'DC3: Step 1 of 9');

fnData = ['data\' buildingName '\dataset.csv'] ;
fnNodes = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesign, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);
seismicCasesIdx = find(ismember(cases, seismicCases));
loading = waitbar(1,loading,'Reading data','Name', 'DC3: Step 1 of 9'); pause(.5);
%%
close(loading); loading = waitbar(0,'Initializing beams','Name', 'DC3: Step 2 of 9'); pause(1);
beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
    
    mAux = [];
    for j = 1 : length(cases)
        M_Ed = DataDesign(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC3beamDesign(fck, fyk , cover, M_Ed, 0);
        mAux(j,:) = [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition];
    end
    [M_rd, conIndex] = max(mAux(:,6)); %best M_rd
    
    sAux = []; sAuxMid = [];
    for j = 1 : length(cases)
        Fz_Ed = DataDesign(barIndex, 4, j);
        given_h = mAux(conIndex, 1);
        given_b = mAux(conIndex, 2);
        longRebarN = mAux(conIndex, 3);
        longRebarPh = mAux(conIndex, 4);
        [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it, sCondition] = DC3beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
        sAux(j,:) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it, sCondition, Fz_Ed];
        
        %midStirrups
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, V_RdMid] = DC3beamDesignMidShear(fck, fyk , cover, Fz_Ed, sec_b, sec_h, longReinfNo);
        sAuxMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, V_RdMid, Fz_Ed];
    end
    
    [V_Rd, conIndex2] = max(sAux(:,4)); %best F_rd
    [V_RdMid, conIndex3] = max(sAuxMid(:,4));
    
    shearPhi = sAux(conIndex2, 1);
    shearSpac = sAux(conIndex2, 2);
    shearLegs = sAux(conIndex2, 3);
    
    shearPhiMid = sAuxMid(conIndex3, 1);
    shearSpacMid = sAuxMid(conIndex3, 2);
    shearLegsMid = sAuxMid(conIndex3, 3);
    
    beams(end+1,:) = [DataDesign(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, mAux(conIndex, 7), shearPhi, shearSpac, shearLegs, V_Rd, sAux(conIndex2, 5), sAux(conIndex2, 6)];
    beamsMid(end+1,:) = [DataDesign(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid, Fz_Ed];
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'DC3: Step 2 of 9');
end
save([folder '\DC3beamsIt1.mat'],'beams')
save([folder '\DC3beamsIt1mid.mat'],'beamsMid')
%%
close(loading); loading = waitbar(0,'Initializing columns','Name', 'DC3: Step 3 of 9'); pause(1);
noStories = max(stories(:,1)); count = 0; columns = []; columnsMid = [];
for i = 1 : size(barsOfColumns,1)
    % #1 design individually bars of a column
    barNames = []; % to append in the end!!!
    for j = 1 : noStories %design of all bars of a column!
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesign(:,1,1) == barName);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        for k = 1 : length(cases)
            N_axial = DataDesign(barIndex, 2, k);
            My_h = DataDesign(barIndex, 5, k);
            Mz_b = DataDesign(barIndex, 6, k);
            V_Ed = DataDesign(barIndex, 4, k);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1(k,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
        end
        %best iteration for that bar
        [~, index] = max(mAux1(:,1));
        %best individual bars of that column
        bestIndi(j,:) = mAux1(index,:);
    end
    
    % #2 design from bottom to top - based on the biggest width of the best individuals
    %bottom bar
    
    bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration
    barName = barsOfColumns(i,1);
    barIndex = find(DataDesign(:,1,1) == barName);
    for p = 1 : length(cases)
        N_axial = DataDesign(barIndex, 2, p);
        My_h = DataDesign(barIndex, 5, p);
        Mz_b = DataDesign(barIndex, 6, p);
        V_Ed = DataDesign(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth);
        mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
    end
    [~, index] = max(mAux3(:,6));%best iteration for that bar
    auxColumns = mAux3(index,:);
    %other bars
    try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
    gWidth = max(minWidth, min(floor(mAux3(index,1) * .9 * 20)/20, mAux3(index,1) - .05));
    
    for j = 2 : noStories
        barName = barsOfColumns(i,j);
        barIndex = find(DataDesign(:,1,1) == barName);
        for p = 1 : length(cases)
            N_axial = DataDesign(barIndex, 2, p);
            My_h = DataDesign(barIndex, 5, p);
            Mz_b = DataDesign(barIndex, 6, p);
            V_Ed = DataDesign(barIndex, 4, p);
%             givenLong = columnComp(mAux3(index, 5),'EC8');%%%%%%%%%%%%%%%%%%%
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth);%, givenLong);
            mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
        end
        [~, index] = max(mAux3(:,6));%best iteration for that bar
        auxColumns(j,:) = mAux3(index,:); %appending
        
        try
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, min(floor(mAux3(index,1) * .9 * 20)/20, mAux3(index,1) - .05));
%             givenLong = mAux3(index, 6);
        end
    end
    
    % check if dimensions agree (no bigger dimension on top of the column)
    % if it fails will try another method to design it
    sorted = issorted(auxColumns(:,1),'descend');
    if ~sorted
        disp('First design method failed')
        %         it will take a lot of time but will fix possibly it
        %         this method will design the column so many times as individual
        %         bars of that column and then, per page of the matrix will put a
        %         column designed based on each bar
        for k = 1 : noStories
            width = bestIndi(k,1); % parent bar
            givenLong = columnComp(bestIndi(k,5), 'EC8');
            
            for j = 1 : noStories
                barName = barsOfColumns(i,j);
                barIndex = find(DataDesign(:,1,1) == barName);
                try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
                Gwidth = max(width, minWidth);
                for p = 1 : length(cases)
                    N_axial = DataDesign(barIndex, 2, p);
                    My_h = DataDesign(barIndex, 5, p);
                    Mz_b = DataDesign(barIndex, 6, p);
                    V_Ed = DataDesign(barIndex, 4, p);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, Gwidth);%, givenLong);
                    ratio = bestIndi(j,6) / reinfPercFin;
                    mAux2(p,:) = [ratio, sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
                end
                [~, index] = max(mAux2(:,7));%best iteration for that bar
                bestIndiMethod(j,:,k) = mAux2(index,:);%best individual bars of that column with a seed provided by the parent bar
            end
        end
        
        possDesigns = [];
        for z = 1 : (noStories) % check if superior dimesions are the same or smaller, if is it, add it to possDesigns
            aCheck = issorted(bestIndiMethod(:,2,z),'descend');
            if aCheck == 1
                possDesigns = [possDesigns, z];
            end
        end
        
        if isempty(possDesigns)
            disp(['Cannot design this bars/column: ' num2str(barNames')]);
            [row, col] = size(auxColumns);
            auxColumns = zeros(row, col);
        else
            validDesigns = bestIndiMethod(:,:,possDesigns);
            for z = 1 : size(validDesigns,3)
                eval(z) = mean(bestIndiMethod(:,1,z));
            end
            [~, bestDesign]= min(eval);
            auxColumns = bestIndiMethod(:,[2:end],bestDesign);
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
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'DC3: Step 3 of 9');
end
save([folder '\DC3columnsIt1.mat'],'columns')
save([folder '\DC3columnsIt1mid.mat'],'columnsMid')
%% 1.3 comparison
% updating the MRD values
close(loading); loading = waitbar(0,'Initializing columns resisting bending moments update','Name', 'DC3: Step 4 of 9'); pause(1);
for i = 1 : size(columns,1)
    b = columns(i, 3);
    h = columns(i, 2);
    areaRebar = columns(i, 6);
    
    mAux = [];
    for j = 1 : length(seismicCases)
        N_Axial = DataDesign(DataDesign(:,1,1) == columns(i,1) , 2, seismicCasesIdx(j));
        M_Rd = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
        mAux(end+1) = M_Rd;
    end
    
    columns(i, 8) = min(mAux);%M_Rd;
    loading = waitbar(i / size(columns,1),loading,'Columns resisting bending moments update progress','Name', 'DC3: Step 4 of 9');
end
%%
close(loading); loading = waitbar(0,'Initializing bending comparisons','Name', 'DC3: Step 5 of 9'); pause(1);
increNeed = [];
for i = 1 : size(nodes,1)
    %skip base points and on the top floor the verification is not needed
    if nodes(i, 5) == 0 || nodes(i, 5) == noStories
        continue 
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
        bendRdZ = bendRdZ + columns(columnRow, 8);
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
    loading = waitbar(i / (size(nodes,1)),loading,'Comparisons progress','Name', 'DC3: Step 5 of 9');
end

close(loading); loading = waitbar(0,'Initializing column re-reinforcement','Name', 'DC3: Step 6 of 9'); pause(1);
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
    
    % iterations
    while ~all(table(:,8) < 0)%while there is a "needed" bigger than what it is
        [~, index] = max(table(:, 8));
        improve = table(index, 7) / 2;
        if max(table(index, 3), table(index, 5)) > improve * 1.2 % if one if more than 60% of total needed just improve the weak one
            [val, col] = min([table(index, 3), table(index, 5)]);
            bar = table(index, col * 2);
            mImprove = table(index, 7) - val;
            toGive = [bar mImprove];
        else
            toGive(:,1) = [table(index, 2); table(index,4)];
            toGive(:,2) = [improve; improve];
        end
        % para a mesma taxa de armadura e sendo quadrados os pilares, a
        % situação sísmica condicionante para o somatório é quando o
        % esforço axial é máximo (das combinações sísmicas) pois reduz o
        % momendo resistente
        % se eu dimensionar para o máximo axial e para o momento
        % necessário, só numa direção (hack: relação dos momentos tem de
        for j = 1 : size(toGive,1)
            N_axial = max(DataDesign(DataDesign(:,1,1) == toGive(j,1), 2, seismicCasesIdx));
            My_h = toGive(j, 2); Mz_b = 0;
            [columnsRow,~] = find(columns(:,1) == toGive(j,1));
            gWidth = columns(columnsRow, 2);
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth);
            V_Ed = columns(columnsRow, 15);
            finalColumn = [toGive(j, 1), sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
            
            newColumns = [newColumns; finalColumn];
            
            %update my table!
            [updateRow, updateCol] = find(table(:,[2,4]) == toGive(j,1));
            updateCol = updateCol * 2;
            
            for k = 1 : length(updateCol)
                table(updateRow(k), [updateCol(k), updateCol(k) + 1]) = [toGive(j, 1), M_Rd];
            end
            
            table(:, 6) = table(:, 3) + table(:, 5);
            table(:, 8) = table(:, 7) - table(:, 6);
        end
    end
    
    %add the ones without changes and clean extra bars
    nonReReinf = setxor(aux1,newColumns(:,1));
    for j = 1 : length(nonReReinf)
        mAux = columns(columns(:,1) == nonReReinf(j),:);
        newColumns = [newColumns; mAux];
    end
    
    for j = 1 : noStories
        lastIndex = find(aux1(j) == newColumns(:,1), 1, 'last');
        primeColumns(j,:) = newColumns(lastIndex,:);
    end
    
    clear auxElement aux1 table toGive
    
    %give concrete dimensions
    for z = 1 : (noStories) % check if superior dimesions are the same or smaller, if is it, add it to possDesigns
        aCheck = issorted(primeColumns(:,2),'descend');
        if ~aCheck
            for j = noStories : -1 : 2
                if primeColumns(j-1,2) < primeColumns(j,2)
                    primeColumns(j-1,[2,3]) = primeColumns(j,[2,3]);
                    b = primeColumns(j-1,2); h = b;
                    areaRebar = primeColumns(j-1, 6);
                    N_Axial = max(DataDesign(DataDesign(:,1,1) == primeColumns(j-1,1) , 2, seismicCasesIdx));
                    primeColumns(j-1,7) = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
                end
            end
        end
    end
    
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
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);%, givenLong);
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
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);%, givenLong);
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
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);%, givenLong);
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
                    [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea);%, givenLong);
                    auxRefs(k, [4:6], j) = [phiRebar, noRebar, areaRebar];
                    givenLong = columnComp(noRebar, phiRebar, 'EC8', 'yes');
                end
        end
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
        N_axial = max(DataDesign(DataDesign(:,1,1) == primeColumns(j,1), 2, seismicCasesIdx));
        primeColumns(j, 8) = MrdColumn(fck, fyk, sec_b, sec_h, areaRebar, N_axial);
        
        longReinfN = primeColumns(j, 4);
        longReinfPh = primeColumns(j, 5);
        Ved = columns(columns(:,1) == primeColumns(j, 1), 15);
        [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, Vrd, sCondition] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, sec_b, [longReinfPh,longReinfN,areaRebar], [], sec_b, longReinfN, longReinfPh);
        primeColumns(j, [9:15]) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, Vrd, sCondition, Ved];
    end
    
    columns13 = [columns13; primeColumns];
    
    clear primeColumns; newColumns = [];
    loading = waitbar(i / (size(barsOfColumns,1)),loading,'Column re-reinforcement progress','Name', 'DC3: Step 6 of 9');
end
save([folder '\DC3columnsIt2rule13.mat'],'columns13')
%% seismic equilibrium
%beam
close(loading); loading = waitbar(0,'Initializing Beam seismic equilibrium','Name', 'DC3: Step 7 of 9'); pause(1);
seismicBeams = []; seismicBeamsMidShear = [];
for i = 1 : length(beamDesiOrd)
    % get beam, get floor, get length, get direction
    % for each one get node, for each node get the sum
    barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
    
    MRb = beams(beams(:,1) == beamDesiOrd(i), 6);
    res = element(element(:,1) == beamDesiOrd(i), [2, 3, 4]);
    node1 = res(1); node2 = res(2); direction = res(3);
    auxNode = [node1, node2];
    cLength = DataDesign(barIndex, 7, 1);
    
    Mid = [];
    for k = 1 : 2
        %get bars on that point
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == direction); barsH = element(row,1);
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
        
        %get design bending moments on those bars
        bendRdH = 0;
        for j = 1 : length(barsH)
            [beamRow, ~] = find(beams(:,1) == barsH(j));
            bendRdH = bendRdH + beams(beamRow, 6);
        end
        
        bendRdZ = 0;
        for j = 1 : length(barsZ)
            [columnRow, ~] = find(columns13(:,1) == barsZ(j));
            bendRdZ = bendRdZ + columns13(columnRow, 8);
        end
        
        if bendRdH < bendRdZ
            Mid(k) = G_RD * MRb;
        else
            Mid(k) = G_RD * MRb * (bendRdZ / bendRdH);
        end
    end
    
    %equilibrium
    Vseismic = max(DataDesign(barIndex, 4, seismicCasesIdx));
    Vshear = Vseismic + sum(Mid) / cLength;
    
    % shear reinforcement
    barIndex = find(beams(:,1) == beamDesiOrd(i));
    given_b = beams(barIndex, 3);
    given_h = beams(barIndex, 2);
    longReinfN =  beams(barIndex, 4);
    longReinfPh =  beams(barIndex, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition] = DC3beamDesign(fck, fyk , cover, 10, Vshear, given_b, given_h, longReinfN, longReinfPh);
    seismicBeams = [seismicBeams; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition, Vshear]];
    
    %mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = DC3beamDesignMidShear(fck, fyk , cover, Vshear, given_b, given_h, longReinfN);
    seismicBeamsMidShear = [seismicBeamsMidShear; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, 0, Vshear]];
    
    loading = waitbar(i / (length(beamDesiOrd)),loading,'Beam seismic equilibrium','Name', 'DC3: Step 7 of 9');
end
save([folder '\DC3beamsIt2SeismicEq.mat'],'seismicBeams')
save([folder '\DC3beamsIt2SeismicEqMid.mat'],'seismicBeamsMidShear')

%column
close(loading); loading = waitbar(0,'Initializing column seismic equilibrium','Name', 'DC3: Step 8 of 9'); pause(1);
pilars = setxor(beamDesiOrd, element(:,1));
seismicColumns = []; seismicColumnsMidShear = [];
for i = 1 : length(pilars)
    barRow = find(columns13(:,1) == pilars(i));
    %MRc = columns(barRow, 8);
    res = element(element(:,1) == pilars(i), [2 3]);
    node1 = res(1); node2 = res(2);
    auxNode = [node1, node2];
    cLength = DataDesign(DataDesign(:,1,1) == pilars(i) , 7, 1);
    % para cada nó!! ir calcular o máximo MRD e somá-los, o máximo porque    % também fz com que  agrave, implica o mínimo axial das sísmicas
    % para cada nó, para cada direção ir buscar o máximo porque quanto    % maior for o somatorio nas beams, maior será o momento midbuildingName
    % fazer a conta do Mid    % design to shear
    for k = 1 : 2
        %obtain max sum of bending moment on top each node of the column
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
        bendRdZaux = 0;
        for j = 1 : length(barsZ)
            b = columns13(columns13(:,1) == barsZ(j), 3);
            h = columns13(columns13(:,1) == barsZ(j), 2);
            areaRebar = columns13(columns13(:,1) == barsZ(j), 6);
            N_Axial = min(DataDesign(DataDesign(:,1,1) == barsZ(j) , 2, seismicCasesIdx));
            M_Rd = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
            bendRdZaux = bendRdZaux + M_Rd;
        end
        bendRdZ(k) = bendRdZaux;
        
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
    
    [~, maxIndex] = max([sum(bendRdX), sum(bendRdY)]);
    bendMomsaux = [bendRdX; bendRdY];
    bendMomsBeam = bendMomsaux(maxIndex,:);
    
    b = columns13(barRow, 3);
    h = columns13(barRow, 2);
    areaRebar = columns13(barRow, 6);
    N_Axial = min(DataDesign(DataDesign(:,1,1) == pilars(i) , 2, seismicCasesIdx));
    MRc = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
    
    for j = 1 : 2
        if bendMomsBeam(j) > bendRdZ(j)
            Mid(k) = G_RD * MRc;
        else
            Mid(k) = G_RD * MRc * (bendMomsBeam(j) / bendRdZ(j));
        end
    end
    
    Vshear = sum(Mid) / cLength;
    
    given_h = columns13(barRow, 2);
    longReinfN = columns13(barRow, 4);
    longReinfPh = columns13(barRow, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, given_h, [longReinfPh,longReinfN,areaRebar], Vshear, given_h, longReinfN, longReinfPh);
    seismicColumns = [seismicColumns; [columns13(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, Vshear]];
    
    %mid shear
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = DC3columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar);
    seismicColumnsMidShear = [seismicColumnsMidShear; [columns13(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, 0, 0]];
    
    loading = waitbar(i / (length(pilars)),loading,'Column seismic equilibrium progress','Name', 'DC3: Step 8 of 9');
end
save([folder '\DC3columnsIt3SeismicEq.mat'],'seismicColumns')
save([folder '\DC3columnsIt3SeismicEqMid.mat'],'seismicColumnsMidShear')
%%
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'DC3: Step 9 of 9'); pause(1);
toSeismo(seismicColumns(:,[1:3 , 5, 4, 9:11]), seismicBeams(:,[1:5, 8:10]), nodes, element, stories, fck, fyk, cover, folder);
loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'DC3: Step 9 of 9');
%%
time = toc;
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);