%% DC2frameDesigner
function [] = DC2frameDesigner(buildingName, fck, fyk, cover, seismicCases)


% clear
% clc
% buildingName = 'regular_DC2' ;
% fck = 30 ;
% fyk = 400;
% cover = .035 ;
% seismicCases = [24:31];

%%
G_RD = 1.1; %factor accounting for overstrength due to steel strain hardening and confinement of the concrete of the compression zone of the section
%%
tic; disp('Started');
loading = waitbar(0,'Reading data','Name', 'Step 1 of 6');

fnData = ['data\' buildingName '\dataset.csv'] ;
fnNodes = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesign, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);
seismicCasesIdx = find(ismember(cases, seismicCases));
loading = waitbar(1,loading,'Reading data','Name', 'Step 1 of 8'); pause(.5);
%%
close(loading); loading = waitbar(0,'Initializing beams','Name', 'Step 2 of 8'); pause(.5);
beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
    
    mAux = [];
    for j = 1 : length(cases)
        M_Ed = DataDesign(barIndex, 5, j);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, roMinCondition] = DC2beamDesign(fck, fyk , cover, M_Ed, 0);
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
        [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it, sCondition] = DC2beamDesign(fck, fyk , cover, M_Ed, Fz_Ed, given_b, given_h, longRebarN, longRebarPh);
        sAux(j,:) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd_it, sCondition];
        
        %midStirrups
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, V_RdMid] = DC2beamDesignMidShear(fck, fyk , cover, Fz_Ed, sec_b, sec_h, longReinfNo);
        sAuxMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, V_RdMid];
    end
    
    [V_Rd, conIndex2] = max(sAux(:,4)); %best F_rd
    [V_RdMid, conIndex3] = max(sAuxMid(:,4));
    
    shearPhi = sAux(conIndex2, 1);
    shearSpac = sAux(conIndex2, 2);
    shearLegs = sAux(conIndex2, 3);
    
    shearPhiMid = sAuxMid(conIndex3, 1);
    shearSpacMid = sAuxMid(conIndex3, 2);
    shearLegsMid = sAuxMid(conIndex3, 3);
    
    beams(end+1,:) = [DataDesign(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, mAux(conIndex, 7), shearPhi, shearSpac, shearLegs, V_Rd, sAux(conIndex2, 5)];
    beamsMid(end+1,:) = [DataDesign(barIndex,1,1), given_h, given_b, longRebarN, longRebarPh, M_rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid];
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'Step 2 of 8');
end
save('toSave\beamsIt1.mat','beams')
save('toSave\beamsIt1mid.mat','beamMid')
%%
close(loading); loading = waitbar(0,'Initializing columns','Name', 'Step 3 of 8'); pause(1);
noStories = max(stories(:,1)); count = 0; columns = []; columnsMid = [];
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
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, minWidth);
            mAux1(k,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition];
        end
        [~, index] = max(mAux1(:,6));%best iteration for that bar
        bestIndi(j,:) = mAux1(index,:);%best individual bars of that column
    end
    
    % 2 design from bottom to top - based on the biggest width of the best individuals
    %bottom bar
    bigOrigWidth = max(bestIndi(:,1)); %biggest width on the individual iteration
    barName = barsOfColumns(i,1);
    barIndex = find(DataDesign(:,1,1) == barName);
    for p = 1 : length(cases)
        N_axial = DataDesign(barIndex, 2, p);
        My_h = DataDesign(barIndex, 3, p);
        Mz_b = DataDesign(barIndex, 4, p);
        try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
        gWidth = max(minWidth, bigOrigWidth);
        [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth);
        mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition];
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
            My_h = DataDesign(barIndex, 3, p);
            Mz_b = DataDesign(barIndex, 4, p);
            givenLong = columnComp(mAux3(index, 5),'EC8');
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, givenLong);
            mAux3(p,:) = [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition];
        end
        [~, index] = max(mAux3(:,6));%best iteration for that bar
        auxColumns(j,:) = mAux3(index,:); %appending
        
        try
            try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
            gWidth = max(minWidth, min(floor(mAux3(index,1) * .9 * 20)/20, mAux3(index,1) - .05));
            givenLong = mAux3(index, 6);
        end
    end
    
    % check if dimensions agree (no bigger dimension on top of the column)
    % if it fails will try another method to design it
    sorted = issorted(auxColumns(:,1),'descend');
    if sorted == 0
        disp('First design method failed')
        %         it will take a lot of time but will fix possibly it
        %         this method will design the column so many times as individual
        %         bars of that column and then, per page of the matrix will put a
        %         column designed based on each bar
        for k = 1 : noStories
            width = bestIndi(k,1); %let's call it parent bar
            givenLong = columnComp(bestIndi(k,5), 'EC8');
            
            for j = 1 : noStories
                barName = barsOfColumns(i,j);
                barIndex = find(DataDesign(:,1,1) == barName);
                try [minWidth] = minWidFind(barName, element, beams); catch minWidth = .2; end
                Gwidth = max(width, minWidth);
                for p = 1 : length(cases)
                    N_axial = DataDesign(barIndex, 2, p);
                    My_h = DataDesign(barIndex, 3, p);
                    Mz_b = DataDesign(barIndex, 4, p);
                    [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, Gwidth, givenLong);
                    ratio = bestIndi(j,6) / reinfPercFin;
                    mAux2(p,:) = [ratio, sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition];
                end
                [~, index] = max(mAux2(:,7));%best iteration for that bar
                bestIndiMethod(j,:,k) = mAux2(index,:);%best individual bars of that column with a seed provided by the parent bar
            end
        end
        
        possDesigns = [];
        for z = 1 : (noStories) % check if superior dimesions are the same or smaller, if is it, add it to possDesigns
            aCheck = issorted(bestIndiMethod(:,1,z),'descend');
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
            auxColumns = bestIndiMethod([2:end],:,bestDesign)
        end
    end
    
    for j = 1 : size(auxColumns,1)
        sec_h = auxClumns(j, 1);
        noRebar = auxClumns(j, 3);
        phiRebar = auxClumns(j, 4);
        [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid] = DC2columnDesignMidShear(fck, fyk , cover, sec_h, noRebar, phiRebar);
        auxColumnsMid(j,:) = [shearReinfPhiMid, shearReinfSpacMid, shearReinfLoopsMid, shearReinfAreaMid];
    end
    
    columns = [columns; [barNames, auxColumns]];
    columnsMid = [columnsMid; [barNames, auxColumns(:,[1:7]), auxColumnsMid]];
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'Step 3 of 8');
end
save('toSave\columnsIt1.mat','columns')
save('toSave\columnsIt1mid.mat','columnsMid')
%% 1.3 comparison
close(loading); loading = waitbar(0,'Initializing bending comparisons','Name', 'Step 4 of 6'); pause(1);
increNeed = [];
for i = 1 : size(nodes,1)
    if nodes(i, 5) == 0 || nodes(i, 5) == noStories
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
        bendRdZ = bendRdZ + columns(columnRow, 6);
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
        
        auxBendMomY = 0;
    end
    increNeed(end+1,:) = [nodes(i,1), max(auxBendMomX, auxBendMomY)];
    loading = waitbar(i / (size(nodes,1) * (noStories - 1)),loading,'Comparisons progress','Name', 'Step 4 of 8');
end

close(loading); loading = waitbar(0,'Initializing column re-reinforcement','Name', 'Step 5 of 8'); pause(1);
newColumns = [];
for i = 1 : size(barsOfColumns,1)
    % construir tabela
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
    while all(table(:,8) >= 0)%houver uma diff positiva
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
        % ser dez, desta forma o programa assume só num sentido)
        for j = 1 : size(toGive,1)
            N_axial = max(DataDesign(DataDesign(:,1,1) == toGive(j,1), 2, seismicCasesIdx));
            My_h = toGive(j, 2);
            Mz_b = My_h / 10;
            [columnsRow,~] = find(columns(:,1) == toGive(j,1))
            gWidth = columns(columnsRow, 2);
            input1 = columns(columnsRow, 4);
            input2 = columns(columnsRow, 5);
            g_reinf = columnComp(input1, input2, 'EC8', 'yes');
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd] = DC2columnDesign(fck, fyk , cover, N_axial, My_h, Mz_b, gWidth, g_reinf);
            if sec_h ~= gWidth; disp("After 1.3 re-reinforcement dimensions don't match!"); end
            shearReinfPhi = columns(columnsRow, 9);
            shearReinfSpac = columns(columnsRow, 10);
            shearReinfLoops = columns(columnsRow, 11);
            shearReinfArea = columns(columnsRow, 12);
            finalColumn = [toGive(j, 1), sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
            newColumns = [newColumns; finalColumn];
            %update my table!
            [updateRow, updateCol] = find(table == toGive(j,1));
            for k = 1 : size(updateCol, 1)
                table(updateRow(k), [updateCol(k), updateCol(k) + 1]) = [toGive(j, 1), M_Rd];
            end
        end
    end
    clear aux1 auxElement table
    loading = waitbar(i / (size(barsOfColumns,1)),loading,'Column re-reinforcement progress','Name', 'Step 5 of 8');
end

%update
for i = 1 : size(newColumns, 1)
    columns(columns(:,1) == newColumns(i,1),:) = newColumns(i,:);
end
save('toSave\columnsIt2rule13.mat','columns')
%% seismic equilibrium
close(loading); loading = waitbar(0,'Initializing Beam seismic equilibrium','Name', 'Step 6 of 8'); pause(1);
seismicBeams = [];
% perimeter = [max(nodes(:,2), max(nodes(:,3), min(nodes(:,2), min(nodes(:,3)];
for i = 1 : length(beamDesiOrd)
    % get beam, get floor, get length, get direction
    % for each one get node, for each node get the sum
    barIndex = find(DataDesign(:,1,1) == beamDesiOrd(i));
    
    MRb = beams(beams(:,1) == beamDesiOrd(i), 6);
    [node1, node2, direction]= element(element(:,1) == beamDesiOrd(i), [2 3 4]);
    auxNode = [node1, node2];
    %     [x1, y1] = nodes(nodes(:,1) == node1,[2, 3]);
    %     [x2, y2, floor] = nodes(nodes(:,1) == node2,[2, 3, 5]);
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
            [columnRow, ~] = find(columns(:,1) == barsZ(j));
            bendRdZ = bendRdZ + columns(columnRow, 6);
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
    longRebarN =  beams(barIndex, 4);
    longRebarPh =  beams(barIndex, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition] = DC2beamDesign(fck, fyk , cover, 10, Vshear, given_b, given_h, longRebarN, longRebarPh);
    seismicBeams = [seismicBeams; [beams(barIndex,[1:7]),shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd, sCondition]];
    loading = waitbar(i / (length(beamDesiOrd)),loading,'Beam seismic equilibrium','Name', 'Step 6 of 8');
end
save('toSave\beamsIt2SeismicEq.mat','seismicBeams')

close(loading); loading = waitbar(0,'Initializing column seismic equilibrium','Name', 'Step 7 of 8'); pause(1);
pilars = setxor(beamDesiOrd, element(:,1));
seismicColumns = [];
for i = 1 : length(pilars)
    barRow = find(columns(:,1) == pilars(i));
    MRc = columns(barRow, 8);
    [node1, node2]= element(element(:,1) == pilars(i), [2 3]);
    auxNode = [node1, node2];
    cLength = DataDesign(DataDesign(:,1,1) == pilars(i) , 7, 1);
    
    % para cada nó!! ir calcular o máximo MRD e somá-los, o máximo porque
    % também fz com que  agrave, implica o mínimo axial das sísmicas
    
    % para cada nó, para cada direção ir buscar o máximo porque quanto
    % maior for o somatorio nas beams, maior será o momento midbuildingName
    %
    % fazer a conta do Mid
    % design to shear
    
    for k = 1 : 2
        %obtain max sum of bending moment on top each node of the column
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 3); barsZ = element(row,1);
        bendRdZaux = 0;
        for j = 1 : 2
            b = columns(columns(:,1) == barsZ(j), 3);
            h = columns(columns(:,1) == barsZ(j), 2);
            areaRebar = columns(columns(:,1) == barsZ(j), 6);
            N_Axial = min(DataDesign(DataDesign(:,1,1) == barsZ(j) , 1, seismicCasesIdx));
            M_Rd = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial);
            bendRdZaux = bendRdZaux + M_Rd;
        end
        bendRdZ(k) = bendRdZaux;
        
        %obtain max sum of bending moment on top each node of the beams, each direction
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 1); barsX = element(row,1);
        bendRdXaux = 0;
        for j = 1 : 2
            M_Rd = beams(beams(:,1) == barsX(j), 6);
            bendRdXaux = bendRdXaux + M_Rd;
        end
        bendRdX(k) = bendRdXaux;
        
        [row, ~] = find(element(:,[2 3]) == auxNode(k) & element(:,4) == 2); barsY = element(row,1);
        bendRdYaux = 0;
        for j = 1 : 2
            M_Rd = beams(beams(:,1) == barsY(j), 6);
            bendRdYaux = bendRdYaux + M_Rd;
        end
        bendRdY(k) = bendRdYaux;
    end
    
    [~, maxIndex] = max([sum(bendRdX), sum(bendRdY)]);
    bendMomsaux = [bendRdX; bendRdY];
    bendMomsBeam = bendMomsaux(maxIndex,:);
    
    for j = 1 : 2
        if bendMomsBeam(j) > bendRdZ(j)
            Mid(k) = G_RD * MRc;
        else
            Mid(k) = G_RD * MRc * (bendMomsBeam(j) / bendRdZ(j));
        end
    end
    
    Vshear = sum(Mid) / cLength;
    
    given_h = columns(barRow, 2);
    longReinfN = columns(barRow, 4);
    longReinfPh = columns(barRow, 5);
    [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC2columnDesign(fck, fyk , cover, 10, 10, 10, 10, [1:3], given_h, longReinfN, longReinfPh)
    seismicColumns = [seismicColumns; [columns(barRow,[1:8]),shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition]];
    loading = waitbar(i / (length(pilars)),loading,'Beam seismic equilibrium','Name', 'Step 7 of 8');
end
save('toSave\columnsIt3SeismicEq.mat','seismicColumns')
%%
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'Step 8 of 8'); pause(1);

% loading = waitbar(0.5, loading,'Updating and writing to Seismo', 'Step 1 of ','Name', 'Step 8 of 6');

toSeismo(columns(:,[1:12]), beams(:,[1:5, 7:11]), nodes, element, stories, fck, fyk, cover)

loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'Step 8 of 8');

time = toc;
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);