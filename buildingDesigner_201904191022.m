%Building Design
clear
clc
tic
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
row = 1;
rowz = 1;
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


%%cagar na eficiência
%fazer acontecer

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
end

%%get the load envelope for each BEAM
row = 1;
rowz = 1;
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

%precedencia de beams para aqueles que intersetam middle way





%{


%bars of the same column %% zero efficiency writes the same line as the
%same number of stories

% problem case: column missing in some stories are not the same column AKA
% not same line
´
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

%}



% no caso de ser a mesma viga, substituir o nome da barra pela de menos
% índice e guardar essa informação
%
%for each element for each node check elements that converge there for each
%check its direction, if same direction no vertical or verticl element but
%its upwardds assign same name and do the the data sorting again but before
%save which bars are same column
%
%


%
    

%{




reinfSolu = zeros(noBars, 29);

%}
%clear aStory auxMat auxMax auxNode auxStories barNameAux bocAux i ind j k ...
%    noCol node_i node_j norm normVect noRow noTimesNaming row rowz rowzz ...
%    sortData unsData vect xBari xCol xi xj yBari yCol yi yj zBari ZBarj ...
%    zi zj

%%
%design columns iteratively on both directions

%{
%
%% EC 2 Design

for i = 1 : noBars
    
    %longitudinal reinforcement
    b = max(.2, round(sqrt(finalData(i,2) / (fcd * 1000)) * 20)/20);
    diff = -1;
    while diff < 0
        redAxial = finalData(i, 2) / ( b^2 * fcd * 1000);
        benMom = sqrt(finalData(i, 5)^2 + finalData(i, 6)^2);
        redBenMom = benMom / (b^3 * fcd * 1000);
        for j = 1 : 13
            mAux(1,j) = abaco(1,j) - redAxial;
        end
        [value col]= min(abs(mAux));
        for k = 1 : 101
            mAux2(k,1) = abaco(k,1) - redBenMom;
        end
        [value row] = min(abs(mAux2));
        reinfPerc = abaco(row, col);
        Asmin = max([.1 * finalData(i,2) / (fyd * 1000), .002 * b^2]);
        AsAbacus = reinfPerc * b^2 * fcd / fyd;
        reinfArea = max(Asmin, AsAbacus * 2);
    
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfArea > 0
                diffAux(j,1) = longReinforce(j,3) - reinfArea;
            else
                diffAux(j,1) = 1000;
            end
        end
    
        [minDiff, minIndex] = min(diffAux);
        noRebar = longReinforce(minIndex,2);
        phiRebar = longReinforce(minIndex,1);
    
        diff = b - 2 * (cover + .02) - (phiRebar/1000 * (noRebar/4)) - (max([.02, phiRebar/1000]) * (noRebar/4 - 1));
        b = b + .05;
    end
    b = b - .05;
    
    %stirrups
    %phi 8mm 
    
    spacing = min([15 * phiRebar/1000, b, .3]);
    
    for k = 1 : 7
        if sAvaila(k) - spacing <= 0
            auxMatrix(k,1) = sAvaila(k) - spacing;
        else
            auxMatrix(k,1) = -1000;
        end
    end
    
    [value index] = max(auxMatrix);
    
    %final data
    reinfSolu(i,1) = finalData(i,1);
    reinfSolu(i,2) = b;
    reinfSolu(i,3) = noRebar;
    reinfSolu(i,4) = phiRebar;
    
    reinfSolu(i,5) = 8;
    reinfSolu(i,6) = sAvaila(index);
    
    reinfSolu(i,8) = longReinforce(minIndex,3) * (b - 2 * (cover + .02)) * fyd * 1000
    
    %reinfSolu(i,8) = longReinforce(minIndex,3) * fyd / (b * b * fcd);
  %  reinfSolu(i,9) = finalData(i,2) / (b * b * fcd * 1000);
   % reinfSolu(i, 10) = finalData(i,5) / ( b * b * b * fcd * 1000);
   % reinfSolu(i, 11) = finalData(i,6) / ( b * b * b * fcd * 1000);
end
%%

%% EC 2 design

for i = 1 : noBars
    reinfSolu(i,1) = finalData(i,1);
        
    %long rebar
    M_Ed = finalData(i,5);
    h = .25;
    start_h = h;
    ratio = 2;
    while ratio > 1.05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = floor(.6 * h * 20) / 20;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = h - (cover + .02);  %approximated
        
        %long rebar starting values
        redBendMom = M_Ed / (b * d^2 * fcd * 1000);
        reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
        reinfAreaAux = reinfPerc * b * d * fcd / fyd;
        reinfArea = max([reinfAreaAux, .26 * fctm * b * d / fyk, .0013 * b * d]);
        
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfArea > 0
                diffAuxL(j,1) = longReinforce(j,3) - reinfArea;
            else
                diffAuxL(j,1) = 1000;
            end
        end
            [minDiff, minIndexA] = min(diffAuxL);
            noRebar = longReinforce(minIndexA,2);
            phiRebar = longReinforce(minIndexA,1);
        
        
        %long rebar iterations
        clearance = min([3 , (noRebar - 1)]) * .05;
        M_Rd = 0;
        if M_Ed > M_Rd
            for j = 1 : size(longReinforce,1) %doesn't make sense on the first while iteration, but it comes to use from second on
                if longReinforce(j,3) - reinfArea > 0
                    diffAuxL(j,1) = longReinforce(j,3) - reinfArea;
                else
                    diffAuxL(j,1) = 1000;
                end
                [minDiff, minIndex] = min(diffAuxL);
                reinfSolu(i,4) = longReinforce(minIndex,2);
                reinfSolu(i,5) = longReinforce(minIndex,1);
                reinfSolu(i,6) = longReinforce(minIndex,3);
                %check asmax
            end
           
            if reinfSolu(i,4) <= 4
                noRebar = reinfSolu(i,4);
                space = noRebar - 1;
                clearance = noRebar * .05;
            else
                noRebar = reinfSolu(i,4);
                space = ceil(noRebar / 2) - 1;
                clearance = space * .05;
            end
            
            diff = b - 2 * cover - clearance - noRebar * (phiRebar / 1000);
            while diff < 0
                 b = b + .05;
                 diff = b - 2 * cover - clearance - .02 - noRebar * (phiRebar / 1000); %* reinfSol(i,8) 
            end
            
            bOh = b / h;
            while bOh > .8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                h = h + .05;
                bOh = b / h;
            end
            
            d = h - (cover + .02);
            reinfPerc = reinfSolu(i, 6) * fyd / (b * d * fcd);
            redBendMom = interp1(abaco(:,2).', abaco(:,1).', max([.004, reinfPerc]));
            M_Rd = redBendMom * b * d^2 * fcd *1000;
            
        end
        
        final_h = h;
        ratio = final_h / start_h;
        reinfSolu(i,2) = h;
        reinfSolu(i,3) = b;
        h = ceil((final_h * .2 + start_h * .8) * 20) / 20;
        start_h = h;
        
    end
    reinfSolu(i, 11) = reinfSolu(i, 6) * (reinfSolu(i, 2) - 2 * (cover + .02)) * fyd * 1000;
    
    %stirrups
    z = h - 2 * (cover + .02);  %approximated 
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5)); %p.100 %9.4 '+' 9.5 EC2 %assuming cot(theta) = 2.5
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 %& (mod(shearReinforce(k,2),2) == mod(space + 1,2) | mod(shearReinforce(k,2) / 2,2) == mod(space + 1,2))
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    [minDiffS, minIndexS] = min(diffAuxS,[],1);
    reinfSolu(i,7) = shearReinforce(minIndexS, 1);
    reinfSolu(i,8) = shearReinforce(minIndexS, 3);
    reinfSolu(i,9) = shearReinforce(minIndexS, 2);
    reinfSolu(i,10) = shearReinforce(minIndexS, 4);
    
    redAxial = finalData(i, 2) / (reinfSolu(i, 2) * reinfSolu(i, 3) * fcd * 1000);
    if redAxial > .1
        reinfSolu(i,:) = 666; % string('Not a beam');
        reinfSolu(i,1) = finalData(i,1);
    end
    
end


clear Asw_s b bOh clearance d diff diffAuxL diffAuxS final_h h i j k M_Ed ...
    M_Rd minDiff minDiffS minIndex minIndexA minIndexS noRebar phiRebar ...
    ratio redAxial redBendMom reinfArea reinfAreaAux reinfPerc rowz space start_h z abaco
% adicionar estribos para garantir amarração de 150mm
%SB7zr4RPp
%% EC8 design

roMin = .5 * fctm / fyk ;
sAvaila = unique(shearReinforce(:,3));
for i = 1: noBars
    %long rebar
    reinfAreaEC8 = roMin * reinfSolu(i,2) * reinfSolu(i,3);
    if reinfAreaEC8 > reinfSolu(i,6)
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfAreaEC8 > 0
                diffAux(j,1) = longReinforce(j,3) - reinfAreaEC8;
            else
                diffAux(j,1) = 1000;
            end
        end
        [minDiff, minIndexA] = min(diffAux);
        
        if longReinforce(minIndexA,3) > reinfSolu(i, 6)
            reinfSolu(i, 12) = longReinforce(minIndexA, 2);
            reinfSolu(i, 13) = longReinforce(minIndexA, 1);
            reinfSolu(i, 14) = longReinforce(minIndexA, 3);
        else
            reinfSolu(i, 12) = reinfSolu(i, 4);
            reinfSolu(i, 13) = reinfSolu(i, 5);
            reinfSolu(i, 14) = reinfSolu(i, 6);
        end
    else
        reinfSolu(i, 12) = reinfSolu(i, 4);
        reinfSolu(i, 13) = reinfSolu(i, 5);
        reinfSolu(i, 14) = reinfSolu(i, 6);
    end
    reinfSolu(i, 15) = reinfSolu(i, 14)/ reinfSolu(i, 6) - 1; %increase of steel area in %
    
    %stirrups DC3
    sMatrix = [0 reinfSolu(i, 8); 1 reinfSolu(i, 2) / 4; 2 reinfSolu(i, 7) * .024; 3 .008 * reinfSolu(i, 5)];
    [sMin sIndex] = min(sMatrix(:,2));
    reinfSolu(i, 21) = sMatrix(sIndex, 1);
    for h = 1 : 7
        if sAvaila(h) - sMatrix(sIndex,2) <= 0
            anotherMatrix(i,1) = sAvaila(h) - sMatrix(sIndex,2);
        else
            anotherMatrix(i,1) = -1000;
        end
    end
    
    [sMaxReal sIndexReal] = max(anotherMatrix);
    
    z = reinfSolu(i, 2) - 2 * (cover + .02);
    sw = sAvaila(sIndexReal);
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5));
        
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 & shearReinforce(k,3) == sw
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    
    [minDiffS, minIndexS] = min(diffAuxS);
    reinfSolu(i,17) = shearReinforce(minIndexS, 1);
    reinfSolu(i,18) = shearReinforce(minIndexS, 3);
    reinfSolu(i,19) = shearReinforce(minIndexS, 2);
    reinfSolu(i,20) = shearReinforce(minIndexS, 4);
    reinfSolu(i,22) = reinfSolu(i,20) / reinfSolu(i,10) - 1;
    
  
    %stirrups DC2
    sMatrix = [0 reinfSolu(i, 8); 1 reinfSolu(i, 2) / 4; 2 reinfSolu(i, 7) * .03; 3 .012 * reinfSolu(i, 5)];
    [sMin sIndex] = min(sMatrix(:,2));
    reinfSolu(i, 28) = sMatrix(sIndex, 1);
    for h = 1 : 7
        if sAvaila(h) - sMatrix(sIndex,2) <= 0
            anotherMatrix(i,1) = sAvaila(h) - sMatrix(sIndex,2);
        else
            anotherMatrix(i,1) = -1000;
        end
    end
    
    [sMaxReal sIndexReal] = max(anotherMatrix);
    
    z = reinfSolu(i, 2) - 2 * (cover + .02);
    sw = sAvaila(sIndexReal);
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5));
        
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 & shearReinforce(k,3) == sw
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    
    [minDiffS, minIndexS] = min(diffAuxS);
    reinfSolu(i,24) = shearReinforce(minIndexS, 1);
    reinfSolu(i,25) = shearReinforce(minIndexS, 3);
    reinfSolu(i,26) = shearReinforce(minIndexS, 2);
    reinfSolu(i,27) = shearReinforce(minIndexS, 4);
    reinfSolu(i,29) = reinfSolu(i,27) / reinfSolu(i,10) - 1;
end

clear anotherAux auxIndex auxMin cover diffAuxL fcd fck fctm fyd fyk fywd...
    i j k longReinforce minDiff minIndexA noBars noCol noRow reinfAreaEC8...
    roMin sAvaila sIndex sMatrix sMatrix2 sMin shearReinforce diffAux diffAuxS...
    h anotherMatrix minDiffS minIndexS sIndexReal sMaxReal sw z Asw_s


%}
toc