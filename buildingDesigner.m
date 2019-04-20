%Building Design
%clear
%clc
%tic
%%
% Top View
% ^Y(2)
% |
% |
% |
% |
% O - - - - - > 
% Z(3)       X(1)
%% import data and inputs
%input('Filename: ');%%%%%%%%%%%%%%%%%
%data = 'dataset.csv' ; 
unsData = importdata('data\dataset.csv');
element = importdata('data\connectivity.csv');
nodes = importdata('data\nodes.csv');

%abaco = importdata('abacusC12_50S500A1.csv');
%longReinforce = importdata('steel_beam.csv');
%shearReinforce = importdata('steel_shear.csv');

%fck = 25;  %MPa
%fyk = 500;  %MPa
%cover = .035; 

%fcd = fck / 1.5; 
%if fck <= 50
%    fctm = .3 * fck^(2/3);
%else
%    fcm = fck + 8;
%    fctm = 2.12 * log(1 + fcm / 10);
%end

%fyd = fyk / 1.15;
%fywd = fyd; 
%% transform data
%extract only 'real' data
[noRow, noCol] = size(unsData.textdata);
barNameAux = extractBefore(unsData.textdata(2:noRow,:),"/");
for i = 1 : (noRow - 1)
    auxMat(i, 1) = str2double(barNameAux(i, 1));
end
sortData = abs([auxMat, unsData.data]);

%get the load envelope for each bar
row = 1; rowz = 1;
finalData = zeros(size(unique(sortData(:, 1)), 1), 7);
for j = unique(sortData(:, 1)).'
    for i = 1 : (noRow - 1) 
        if j == sortData(i, 1)
        auxBlock(row, :) = sortData(i, :);
        row = row + 1;
        end
    end
    auxMax = max(auxBlock); %auxMax = max(auxBlock,[],1)
    finalData(rowz, :) = auxMax;
    clear auxBlock; % to prevent values of previous bar being there
    rowz = rowz + 1;
    row = 1;
end
noBars = rowz - 1;

%identify and assign stories to nodes
auxStories = unique(nodes(:, 4));
ind = 0;
for k = 1 : size(auxStories)
    stories(ind + 1, 2) = auxStories(k);
    stories(ind + 1, 1) = ind;
    ind = ind + 1;
end

for j = 1 : size(nodes,1)
    nodes(j,5) = stories(stories(:,2) == nodes(j,4), 1);
end

%identify and assign directions to bars
for i = 1 : noBars
    node_i = element(i,2);
    node_j = element(i,3);
    xi = nodes(nodes(:,1) == node_i, 2);
    yi = nodes(nodes(:,1) == node_i, 3);
    zi = nodes(nodes(:,1) == node_i, 4);
    xj = nodes(nodes(:,1) == node_j, 2);
    yj = nodes(nodes(:,1) == node_j, 3);
    zj = nodes(nodes(:,1) == node_j, 4);
    vect = [xi-xj, yi-yj, zi-zj];
    norm = sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);
    normVect = abs(vect / norm);
    if normVect == [1, 0, 0]
        element(i, 4) = 1;
    elseif normVect == [0, 1, 0]
        element(i, 4) = 2;
    elseif normVect == [0, 0, 1]
        element(i, 4) = 3;
    end
end

%bars to beams
barsOfBeams =[];
for j = 1 : size(element,1)
    for i = [2 3]
        auxNode = element(j,i);
        row = 1;
        for k = 1 : size(element,1)
            if element(k, 2) == auxNode | element(k, 3) == auxNode
                auxAux(row,1) = element(k,1);
                auxAux(row,2) = element(k,4);
                auxAux(row,3) = element(k,2);
                auxAux(row,4) = nodes(nodes(:,1) == auxAux(row,3), 5);
                auxAux(row,5) = element(k,3);
                auxAux(row,6) = nodes(nodes(:,1) == auxAux(row,5), 5);
                row = row + 1;
            end
        end
        
        dir1 = sortrows(auxAux(auxAux(:,2) == 1,:));
        dir2 = sortrows(auxAux(auxAux(:,2) == 2,:));
        dir3 = sortrows(auxAux(auxAux(:,2) == 3,:));
        
        if isempty(dir1) == 0
            if (size(dir1,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir1(1,6)))
                barsOfBeams = [barsOfBeams; (dir1(:,1))'] ;
            end
        end
        
        if isempty(dir2) == 0
            if (size(dir2,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir2(1,6)))
                barsOfBeams = [barsOfBeams; (dir2(:,1))'] ;
            end
        end
        
        clear auxAux dir1 dir2 dir3;
	end
end

%beams ON beams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%new block of code   
%precedencia de beams para aqueles que intersetam middle way
beamsOnBeams =[];
for j = 1 : size(element,1)
    for i = [2 3]
        auxNode = element(j,i);
        row = 1;
        for k = 1 : size(element,1)
            if element(k, 2) == auxNode | element(k, 3) == auxNode
                auxAux(row,1) = element(k,1);
                auxAux(row,2) = element(k,4);
                auxAux(row,3) = element(k,2);
                auxAux(row,4) = nodes(nodes(:,1) == auxAux(row,3), 5);
                auxAux(row,5) = element(k,3);
                auxAux(row,6) = nodes(nodes(:,1) == auxAux(row,5), 5);
                row = row + 1;
            end
        end
        
        dir1 = sortrows(auxAux(auxAux(:,2) == 1,:));
        dir2 = sortrows(auxAux(auxAux(:,2) == 2,:));
        dir3 = sortrows(auxAux(auxAux(:,2) == 3,:));
        
        if size(dir1,1) == 1
            if (size(dir2,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir1(1,6)))
                beamsOnBeams = [beamsOnBeams; [(dir1(:,1))' , dir2(:,1)']] ;
            end
        end
        
        if size(dir2,1) == 1
            if (size(dir1,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir2(1,6)))
                beamsOnBeams = [beamsOnBeams; [(dir2(:,1))' , dir1(:,1)']] ;
            end
        end
        
        clear auxAux dir1 dir2 dir3;
	end
end

%naming
barsOfBeams = unique(barsOfBeams,'rows');
noTimesNaming = max(sum([barsOfBeams(:,1);barsOfBeams(:,2)]==[barsOfBeams(:,1);barsOfBeams(:,2)]'));

for i = 1 : noTimesNaming
    for j = 1 : noBars
        for k = 1 : size(barsOfBeams,1)
            if finalData(j,1) == barsOfBeams(k,1)
                finalData(j,1) = barsOfBeams(k,2);
                %finalData(j,7) = finalData(j,7) + finalData(finalData(j,7) == barsOfBeams(k,1),7)
            end
        end
    end
    
    for a = 1 : size(beamsOnBeams,1)
       for s = 1 : size(beamsOnBeams,2)
           for d = 1 : size(barsOfBeams)
               if beamsOnBeams(a,s) == barsOfBeams(d,1)
                   beamsOnBeams(a,s) = barsOfBeams(d,2);
               end
           end
       end
    end
end

beamsOnBeams = unique(beamsOnBeams,'rows') ;
beamsOnBeams = beamsOnBeams.' ;
beamsOnBeams = unique(beamsOnBeams,'rows') ;
beamsOnBeams = beamsOnBeams.' ;

%order of beams to be designed 
nRow = size(beamsOnBeams,1);
nCol = size(beamsOnBeams,2);
nColNew = 3;
while nColNew > nCol
    for k = 1 : (nColNew - 1) %% inneficcient rewrites same line
        for i = 1 : nRow
            for j = 1 : nRow
                if beamsOnBeams(i, k+1) == beamsOnBeams(j,1)
                    beamsOnBeams(i, k+2 ) = beamsOnBeams(j, 2);
                    %beamsOnBeams(j,:) = [];
                end
            end
        end
    end
    nCol = nColNew;
    nColNew = size(beamsOnBeams,2);
end

revOrder = [];
for i = 1 : size(beamsOnBeams)
    revOrder = [revOrder, beamsOnBeams(i,:)] ;
end
revOrder = unique(revOrder,'stable');
revOrder = revOrder(revOrder ~= 0) ;

beamDesiOrd = [];
for i = size(revOrder,2): -1 : 1
    beamDesiOrd = [beamDesiOrd, revOrder(i)] ;
end

%%get the load envelope for each BEAM
row = 1; rowz = 1;
DataDesign = zeros(size(unique(finalData(:, 1)), 1), 7);
for j = unique(finalData(:, 1)).'
    for i = 1 : noBars
        if j == finalData(i, 1)
        auxBlock(row, :) = finalData(i, :);
        row = row + 1;
        end
    end
    auxMax = max(auxBlock,[],1);
    DataDesign(rowz, :) = auxMax;
    clear auxBlock; % to prevent values of previous bar being there
    rowz = rowz + 1;
    row = 1;
end

%bars of the same column %% zero efficiency writes the same line as the
%same number of stories
% PROBLEM CASE: column missing in some stories are not the same column AKA
% not same line
barsOfColumns = [] ;
bocAux = element(element(:,4) == 3,:) ;
rowzz = 1 ;
for i = 1 : size(bocAux,1)
    xCol = nodes(bocAux(i,2) == nodes(:,1),2);
    yCol = nodes(bocAux(i,2) == nodes(:,1),3);
    for j = 1 : size(bocAux,1)
        xBari = nodes(bocAux(j,2) == nodes(:,1),2);
        yBari = nodes(bocAux(j,2) == nodes(:,1),3);
        zBari = nodes(bocAux(j,2) == nodes(:,1),5);
        %xBarj = nodes(bocAux(j,3) == nodes(:,1),2);
        %yBarj = nodes(bocAux(j,3) == nodes(:,1),3);
        zBarj = nodes(bocAux(j,3) == nodes(:,1),5);
        aStory = max([zBari zBarj]);
        if (xBari == xCol & yBari == yCol)
            barsOfColumns(rowzz, aStory) = bocAux(j,1) ;
        end
    end
    rowzz = rowzz + 1 ;
end                             

clear a aStory auxMat auxMax auxNode auxStories barNameAux ... barsOfBeams ...
    barsOfColumns bocAux d finalData i ind j k nCol nColNew noBars ...
    noCol node_i node_j nodes norm normVect noRow nRow revOrder row ...
    rowz rowzz s sortData unsData vect xBari xCol xi xj yBari yCol yi ...
    yj zBari zBarj zi zj
toc