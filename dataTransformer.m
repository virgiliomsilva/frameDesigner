function [barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, ...
    DataDesign, element, noTimesNaming, stories, nodes, cases] = ...
    dataTransformer (fnData, fnElement, fnNodes)
%% TO DELETE
% tic
% clear
% clc
% % %68zh7q
% % % 
% % fnData = 'data\calcada_da_tapada\datasetall.csv' ;
% % fnNodes = 'data\calcada_da_tapada\nodes.csv' ;
% % fnElement = 'data\calcada_da_tapada\connectivity.csv' ;
% 
% fnData = 'data\regular_DC1\dataset.csv' ;
% fnNodes = 'data\regular_DC1\nodes.csv' ;
% fnElement = 'data\regular_DC1\connectivity.csv' ;
%% IMPORT DATA
unsData = importdata(fnData);
nodes = importdata(fnNodes);
element = importdata(fnElement);

clear fnData fnElement fnNodes
%% ORGANISE DATA
allData = []; cases = [];

for i = 1 : length(unsData)
    mixed = unsData(i,1);
    try
        [BAR, CASE, FX, FY, FZ, MY, MZ, LENGTH] = infoExtractor(mixed);
    catch
        continue
    end
    
    cases(end+1) = CASE;
    
    try
        if allData(end, 7, CASE) == 0
            allData(end, [1:7], CASE)= [BAR, FX, FY, FZ, MY, MZ, LENGTH];
        else
            allData(end+1, [1:7], CASE)= [BAR, FX, FY, FZ, MY, MZ, LENGTH];
        end
    catch
        allData(1, [1:7], CASE)= [BAR, FX, FY, FZ, MY, MZ, LENGTH];
    end
end

cases = unique(cases);
noCases = length(cases);
nameBars = unique(allData(:,1,cases(1)));
noBars = length(nameBars);

clear BAR CASE FX FY FZ i LENGTH mixed MY MZ unsData
%% COMBINATION LOADS
finalData = zeros(noBars, 7, noCases);

for i = 1 : noCases
   for j = 1 : noBars
      finalData(j,:,i) = max(allData(allData(:,1,cases(1)) == nameBars(j),:,cases(i)));
   end
end

clear allData i j nameBars
%% IDENTIFY AND ASSIGN STORIES TO NODES
auxStories = unique(nodes(:, 4));

for i = 1 : length(auxStories)
    stories(i, [1 2]) = [i-1, auxStories(i)] ; % level of fotings is zero
end

for j = 1 : size(nodes,1)
    nodes(j,5) = stories(stories(:,2) == nodes(j,4), 1);
end

clear auxStories i j
%% IDENTIFY AND ASSIGN DIRECTIONS TO BARS
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
    normVect = round(abs(vect / norm));
    if normVect == [1, 0, 0]
        element(i, 4) = 1;
    elseif normVect == [0, 1, 0]
        element(i, 4) = 2;
    elseif normVect == [0, 0, 1]
        element(i, 4) = 3;
    end
end

clear i node_i node_j norm normVect vect xi xj yi yj zi zj
%% BARS to BEAMS & BEAMS on BEAMS
% it works together because it checks each bar two time (one bar = two nodes)
barsOfBeams =[]; beamsOnBeams =[];
for j = 1 : size(element,1)
    for i = [2 3]
        auxNode = element(j,i);
        row = 1;
        for k = 1 : size(element,1)
            if element(k, 2) == auxNode | element(k, 3) == auxNode
                auxAux(row,1) = element(k,1);                              %element
                auxAux(row,2) = element(k,4);                              %direction
                auxAux(row,3) = element(k,2);                              %node 1
                auxAux(row,4) = nodes(nodes(:,1) == auxAux(row,3), 5);     %floor node 1
                auxAux(row,5) = element(k,3);                              %node 2
                auxAux(row,6) = nodes(nodes(:,1) == auxAux(row,5), 5);     %floor node 2
                row = row + 1;
            end
        end
        
        dir1 = sortrows(auxAux(auxAux(:,2) == 1,:));
        dir2 = sortrows(auxAux(auxAux(:,2) == 2,:));
        dir3 = sortrows(auxAux(auxAux(:,2) == 3,:));
        
        if ~isempty(dir1)
            if (size(dir1,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir1(1,6))) % if there is more than one bar in direction 1 AND none in direction 3 OR is a column supported by this beam
                barsOfBeams = [barsOfBeams; (dir1(:,1))'] ;                                             % these bars are the same beam!
            elseif size(dir1,1) == 1
                if (size(dir2,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir1(1,6))) % if there's only one beam in that direction and in one point is supported on another beam
                    beamsOnBeams = [beamsOnBeams; [(dir1(:,1))' , dir2(:,1)']] ;                            % it is supported by that beam
                end
            end
        end
        
        if ~isempty(dir2)
            if (size(dir2,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir2(1,6)))
                barsOfBeams = [barsOfBeams; (dir2(:,1))'] ;
            elseif size(dir2,1) == 1
                if (size(dir1,1) > 1 & (isempty(dir3) == 1 | min(min([dir3(:,4) dir3(:,6)])) >= dir2(1,6)))
                    beamsOnBeams = [beamsOnBeams; [(dir2(:,1))' , dir1(:,1)']] ;
                end
            end
        end
        
        clear auxAux dir1 dir2 dir3;
	end
end

barsOfBeams = unique(barsOfBeams,'rows');

clear auxNode i j k row
%% NAMING - reduce bars to beams
if isempty(barsOfBeams)
    noTimesNaming = 0;
else
    noTimesNaming = max(sum([barsOfBeams(:,1);barsOfBeams(:,2)]==[barsOfBeams(:,1);barsOfBeams(:,2)]'));
    for h = 1 : noCases  
        for i = 1 : noTimesNaming
            for j = 1 : noBars
                for k = 1 : size(barsOfBeams,1)
                    if finalData(j,1,h) == barsOfBeams(k,1)
                        finalData(j,1,h) = barsOfBeams(k,2);
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
    end
end

beamsOnBeams = unique(beamsOnBeams,'rows') ;
beamsOnBeams = beamsOnBeams.' ;
beamsOnBeams = unique(beamsOnBeams,'rows') ; 
beamsOnBeams = beamsOnBeams.' ;

clear a d h i j k s
%% ORDER OF BEAMS TO BE DESIGNED 
beamDesiOrd = [];
noRows = size(beamsOnBeams,1);
if noRows > 0
    nCol = size(beamsOnBeams,2);
    nColNew = nCol + 1;
    while nColNew > nCol
        for k = 1 : (nColNew - 1)
            for i = 1 : noRows
                for j = 1 : noRows
                    
                    try
                        if beamsOnBeams(i, k+1) == beamsOnBeams(j,1)
                            beamsOnBeams(i, k+2 ) = beamsOnBeams(j, 2);
                        end
                    catch
                        continue
                    end
                    
                end
            end
        end
        nCol = nColNew;
        nColNew = size(beamsOnBeams,2);
    end

%     revOrder = [];
%     for i = 1 : size(beamsOnBeams,1)
%         revOrder = [revOrder, beamsOnBeams(i,:)] ;
%     end
%     
%     revOrder = unique(revOrder,'stable'); revOrder = revOrder(revOrder ~= 0) ;
    
%     for i = size(revOrder,2): -1 : 1
%         beamDesiOrd = [beamDesiOrd, revOrder(i)] ;
%     end

    for i = 1 : size(beamsOnBeams,1)
        beamDesiOrd = [beamDesiOrd, fliplr(beamsOnBeams(i,:))];
    end
    
    beamDesiOrd = unique(beamDesiOrd,'stable'); 
    beamDesiOrd = beamDesiOrd(beamDesiOrd ~= 0);
    
    anoAuxBlock = element(element(:,4) ~= 3, 1);
    AA = setxor(beamDesiOrd, anoAuxBlock);
    if isempty(barsOfBeams)
        barsToAdd = AA;
    else
        barsToAdd = setxor(AA,barsOfBeams(:,1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    beamDesiOrd = [beamDesiOrd, barsToAdd'];

else
    
%     anoAuxBlock = element(element(:,4) ~= 3, 1);
%     if isempty(barsOfBeams)
%         barsToAdd = anoAuxBlock;
%     else
%         barsToAdd = setxor(anoAuxBlock,barsOfBeams(:,1)); %%%%%%%%%%%%%%%
%     end


    anoAuxBlock = element(element(:,4) ~= 3, 1);
    if isempty(barsOfBeams)
        beamDesiOrd = anoAuxBlock;
    else
        beamDesiOrd = setxor(anoAuxBlock,barsOfBeams(:,1)); %%%%%%%%%%%%%%%
    end

end

%beamDesiOrd = [beamDesiOrd, barsToAdd'];

 
% tama = size(beamsOnBeams,1);
% 
% for y = 1 : size(barsToAdd,1)
%     beamsOnBeams(y+tama,1) = barsToAdd(y,1) ;
% end

clear AA anoAuxBlock barsToAdd i j k nCol nColNew noRows 
%% LOADS FOR COLUMNS AND FAKE BEAMS
% get the load envelope for each BEAM and COLUMN
DataDesign = [];
for k = 1: noCases
    rowz = 1;
    for j = unique(finalData(:, 1, 1)).'
        row = 1;
        for i = 1 : noBars
            if j == finalData(i, 1, k)
            auxBlock(row, :) = finalData(i, :, k);
            row = row + 1;
            end
        end
        auxMax = max(auxBlock,[],1);
        DataDesign(rowz, :, k) = auxMax;
        clear auxBlock; % prevent values of previous bar being there
        rowz = rowz + 1;
    end
end

clear auxMax noBars noCases finalData i j k row rowz
%% FAKE BEAMS WITH THE SAME NAME AS IN beamDesiOrd EXTREMETY NODES AND REAL BEAM LENGTH
% beams: first and last node and length
auxVec1 = [];
auxVec2 = [];
auxVec3 = [];
fakeBeams = [];

for i = 1 : size(beamDesiOrd,2) % talvz  acrscentar as barras dos pilares
    auxVec1 = [auxVec1, beamDesiOrd(1,i)];
    
    for j = 1 : noTimesNaming %bars of i beam
        if j <= size(auxVec1,2)
        iniBar = barsOfBeams(barsOfBeams(:,2) == auxVec1(1,j),1);
        auxVec1 = [auxVec1; iniBar];
        end
    end
    
    for k = 1 : size(auxVec1,2) % nodes of those bars
       nodz = element(element(:,1) == auxVec1(1,k),[2 3]);
       auxVec2 = [auxVec2, nodz] ;
    end
    auxVec2 = auxVec2' ;
    auxVec2 = unique(auxVec2, 'rows') ;
    auxVec2 = auxVec2' ;
    
    for v = 1 : size(auxVec2,2)
        auxVec3 = [auxVec3; nodes(nodes(:,1) == auxVec2(1,v),[1:4])];
    end
    
    [maxX maxXind]= max(auxVec3(:,2));
    [maxY maxYind]= max(auxVec3(:,3));
  %  [maxZ maxZind]= max(auxVec3(:,4));
    [minX minXind]= min(auxVec3(:,2));
    [minY minYind]= min(auxVec3(:,3));
 %   [minZ minZind]= min(auxVec3(:,4));
    if maxX == minX
        beamLen = maxY - minY;
        nodeI = auxVec3(minYind, 1);
        nodeJ = auxVec3(maxYind, 1);
    elseif maxY == minY
        beamLen = maxX - minX;
        nodeI = auxVec3(minXind, 1);
        nodeJ = auxVec3(maxXind, 1);
    end
    
    fakeBeams = [fakeBeams; [beamDesiOrd(1,i) nodeI nodeJ beamLen]];
    auxVec1 = [];
    auxVec2 = [];
    auxVec3 = [];
end

clear auxVec1 auxVec2 auxVec3 beamLen i iniBar j k maxX maxXind maxY maxYind ...
    minX minXind minY minYind nodeI nodeJ nodz v
%% BARS OF THE SAME COLUMN
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
barsOfColumns = unique(barsOfColumns, 'rows') ;

clear aStory bocAux i j rowzz xBari xCol yBari yCol zBari zBarj
% disp(['** Finished in ', num2str(round(toc,2)), ' secs **']);
% end