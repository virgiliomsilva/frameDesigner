% shotCallerDC3
% beam design long rebar
% store, for each node, for each direction, design bending moment
% design column accordingly max bend moment of each node and robot forces

%% function caller
clear
clc
tic
[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories, nodes] = dataTransformer ();

fck = 25 ;
fyk = 500 ;
fyd = fyk / 1.15;
fcd = fck / 1.5;
fywd = fyd;
cover = .035 ;
hf = .2 ;%slab/flange depth
hfAs =  0.000314159 ;%top reinforcement of the slab per meter in squared meters
gRd = 1.1;
%%
beamsDC3=[];
for i = size(beamsOnBeams,1) : -1 : 1
    for j = size(beamsOnBeams,2) : -1 : 1
        if beamsOnBeams(i,j) ~= 0
            columnMarker = find(beamsOnBeams(i,:)~=0,1,'last');
            if j == columnMarker
                M_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),5) ;
                Fz_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),4) ;
                length = fakeBeams(fakeBeams(:,1) == beamsOnBeams(i,j),4);
                [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd] = beamDesignDC3long(fck, fyk , cover, M_Ed, Fz_Ed, length);
                savedDim = sec_h;
                beamNo = beamsOnBeams(i,j);
                beamsDC3 = [beamsDC3; [beamNo, sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd]];
            else
                M_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),5) ;
                Fz_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),4) ;
                length = fakeBeams(fakeBeams(:,1) == beamsOnBeams(i,j),4);
                [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd] = beamDesignDC3long(fck, fyk , cover, M_Ed, Fz_Ed, length, savedDim);
                beamNo = beamsOnBeams(i,j);
                beamsDC3 = [beamsDC3; [beamNo, sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd]];
            end
        end
    end
    clear columnMarker savedDim 
end
beamsDC3 = unique(beamsDC3,'rows');
beamsDC3(60,:) = [];                %%%%%%%%%%%%%%%%%%%%%%%%

    
bendingNodes = [];
for i = 1 : size(beamsDC3,1)
    node1 = fakeBeams(fakeBeams(:,1) == beamsDC3(i,1),2);
    node2 = fakeBeams(fakeBeams(:,1) == beamsDC3(i,1),3);
    dir = element(element(:,1) == beamsDC3(i,1),4);
    M_Rd = beamsDC3(i,7);
    bendingNodes = [bendingNodes;[node1 dir M_Rd];[node2 dir M_Rd]];
end
%% bending moments accordingly x1.3 rule
aux2 = [];
bending13 = [];
for i = 1 : size(barsOfColumns,1)
    for j = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,j) ~= 0
            nodz1 = element(element(:,1) == barsOfColumns(i,j),[2 3]);
            for k = 1 : 2
                aux = bendingNodes(bendingNodes(:,1) == nodz1(1,k), [2 3]);
                soma = sum(aux(aux(:,1) == 1, 2));
                aux2 = [aux2; soma];
                soma = sum(aux(aux(:,1) == 2, 2));
                aux2 = [aux2; soma];
            end

            newBendMom = max(aux2); % square columns

            if nodes(nodes(:,1) == nodz1(1,1), 5) == max(stories(:,1))
                bending13 = [bending13; [barsOfColumns(i, j), 1.3 * newBendMom]];
            else
                bending13 = [bending13; [barsOfColumns(i, j), 1.3 * newBendMom/2]];
            end
            clear soma 
            aux2 = [];
        end
    end
end
%%
columnsDC3 = [];
bAux = [];
for i = 1 : size(barsOfColumns,1)
    for j = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,j) ~= 0
            N_Axial = DataDesign(DataDesign(:,1) == barsOfColumns(i,j),2);
            My_h = max([bending13(bending13(:,1) == barsOfColumns(i,j),2), DataDesign(DataDesign(:,1) == barsOfColumns(i,j),5)]);
            Mz_b = max([bending13(bending13(:,1) == barsOfColumns(i,j),2), DataDesign(DataDesign(:,1) == barsOfColumns(i,j),6)]);
            [h, b, noRebar, phiRebar] = columnDesignDC3(fck, fyk , cover, N_Axial, My_h, Mz_b);
            bAux = [bAux; [h, b, noRebar, phiRebar]];
        end
    end
    
    hColumn = max(bAux(:,1));
    bColumn = max(bAux(:,2));
    phiColumn = min(bAux(:,4));
    
    for jj = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,jj) ~= 0
            bar = barsOfColumns(i, jj);
            N_Axial = DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),2);
            My_h = max([bending13(bending13(:,1) == barsOfColumns(i,j),2), DataDesign(DataDesign(:,1) == barsOfColumns(i,j),5)]);
            Mz_b = max([bending13(bending13(:,1) == barsOfColumns(i,j),2), DataDesign(DataDesign(:,1) == barsOfColumns(i,j),6)]);
            [h, b, noRebar, phiRebar, areaRebar, M_Rd] = columnDesignDC3(fck, fyk , cover, N_Axial, My_h, Mz_b, bColumn, hColumn, phiColumn);
            columnsDC3 = [columnsDC3; [bar, h, b, noRebar, phiRebar, areaRebar]];
            bendingNodes = [bendingNodes; [element(element(:,1) == barsOfColumns(i,jj),2), 3, M_Rd]; [element(element(:,1) == barsOfColumns(i,jj),3), 3, M_Rd]];
        end
    end
    clear bColumn hColumn phiColumn
    bAux = [];
end
%columns(columns(:,1) == 0,:) = [] ;
%%
%aumentar bending moments às vigas e nós exteriores

abaco = importdata('info\abacusC12_50S500all.mat');

beamsBef = beamDesiOrd.' ; % give a better name to this variable -.-
x_min = min(nodes(:,2));
y_min = min(nodes(:,3));
x_max = max(nodes(:,2));
y_max = max(nodes(:,3));
x_extra = 6.25; 
y_extra = 9.7;

noFlangeBeams = [];
for i = 1 : size(beamsBef,1)
    nodz2 = element(element(:,1) == beamsBef(i,1),[2 3]);
    x_nodei = nodes(nodes(:,1) == nodz2(1,1),2);
    x_nodej = nodes(nodes(:,1) == nodz2(1,2),2);
    y_nodei = nodes(nodes(:,1) == nodz2(1,1),3);
    y_nodej = nodes(nodes(:,1) == nodz2(1,2),3);
    
    if (((x_nodei == x_nodej) & (x_nodei == x_min | x_nodei == x_max | x_nodei == x_extra))) | ...
            (((y_nodei == y_nodej) & (y_nodei == y_min | y_nodei == y_max | y_nodei == y_extra)))
        noFlangeBeams = [noFlangeBeams; beamsBef(i,1)];
        beamsBef(i,1) = 0;
    end
end
beamsBef(beamsBef(:,1) == 0,:) = [];

auxMat = element(element(:,4) == 3, :);
columnWidth1 = 0;
columnWidth2 = 0;
interiorBeams = []; %%%%%%%
for ii = 1 : size(beamsBef,1)
    nodz = fakeBeams(fakeBeams(:,1) == beamsBef(ii,1),[2 3]);
    x_nodei = nodes(nodes(:,1) == nodz(1,1),2);
    x_nodej = nodes(nodes(:,1) == nodz(1,2),2);
    y_nodei = nodes(nodes(:,1) == nodz(1,1),3);
    y_nodej = nodes(nodes(:,1) == nodz(1,2),3);
    
    column1 = auxMat(auxMat(:,2) == nodz(1,1),1); if isempty(column1), column1 = auxMat(auxMat(:,3) == nodz(1,1),1); end; 
    if isempty(column1) == 0, columnWidth1 = columnsDC3(columnsDC3(:,1) == column1(1,1),2); end;
    column2 = auxMat(auxMat(:,2) == nodz(1,2),1); if isempty(column2), column2 = auxMat(auxMat(:,3) == nodz(1,2),1); end; 
    if isempty(column2) == 0, columnWidth2 = columnsDC3(columnsDC3(:,1) == column2(1,1),2); end;
    if (isempty(columnWidth1) & isempty(columnWidth2))
        interiorBeams = [interiorBeams; beamsBef(ii,1)];
        continue;
    end
    minColWid = min([columnWidth1, columnWidth2]);
    
    reinfArea = max(beamsDC3(beamsDC3(:,1) == beamsBef(ii,1),6));
    b = beamsDC3(beamsDC3(:,1) == beamsBef(ii,1),3); %%%%%%%%%%%
    h = beamsDC3(beamsDC3(:,1) == beamsBef(ii,1),2);%%%%%%%%%%%%%
    d = h - (cover + .02);
    
    if x_nodei == x_min | x_nodei == x_max | x_nodei == x_extra | x_nodej == x_min | ... 
            x_nodej == x_max | x_nodej == x_extra | y_nodei == x_min | y_nodei == x_max | ...
            y_nodei == x_extra | y_nodej == x_min | y_nodej == y_max | y_nodej == x_extra
        
        width = minColWid + 2 * hf;
    else
        width = minColWid + 4 * hf;
    end
    
    flangeW = width - beamsDC3(beamsDC3(:,1) == beamsBef(ii,1),3); %%%%%%%%%%%%%%
    ratio = reinfArea / (reinfArea + hfAs * flangeW) ;
    reinfPerc = reinfArea * fyd / ( b * d * fcd);
    redBenMom = abaco(ratio, reinfPerc);
    MRbi = redBenMom * b * d^2 * fcd * 1000;
    %%%%%%%%%%%%%%%%%%%%%%%
    direc = element(element(:,1) == beamsBef(ii,1),4);
    nodz = fakeBeams(fakeBeams(:,1) == beamsBef(ii,1),[2 3]);
    
    auxVec = [];
    for k = 1 : 2
        SmRb = sum(bendingNodes(bendingNodes(:,1) == nodz(1,k) & bendingNodes(:,1) == direc,3));
        SmRc = sum(bendingNodes(bendingNodes(:,1) == nodz(1,k) & bendingNodes(:,1) == 3,3));
        if SmRb > SmRc
            M_Ed = gRd * MRbi * SmRc / SmRc;
        else
            M_Ed = gRd * MRbi;
        end
        auxVec = [auxVec, M_Ed];
    end
    length = fakeBeams(fakeBeams(:,1) == beamsBef(ii,1),4);
    z1 = (auxVec(1) - auxVec(2) + DataDesign(DataDesign(:,1) == beamsBef(ii,1),4)) / length;
    z2 = (auxVec(2) - auxVec(1) + DataDesign(DataDesign(:,1) == beamsBef(ii,1),4)) / length;
    shearForce = max([z1 z2]);
    fakeBeams(fakeBeams(:,1) == beamsBef(ii,1),5) = shearForce ;
    %%%%%%%%%%%%%%%%%%%%%%%
    %clear b h d reinfArea
end

for i = 1 : size(noFlangeBeams,1)
    nodz = fakeBeams(fakeBeams(:,1) == noFlangeBeams(i,1),[2 3]);
    x_nodei = nodes(nodes(:,1) == nodz(1,1),2);
    x_nodej = nodes(nodes(:,1) == nodz(1,2),2);
    y_nodei = nodes(nodes(:,1) == nodz(1,1),3);
    y_nodej = nodes(nodes(:,1) == nodz(1,2),3);
    
    column1 = auxMat(auxMat(:,2) == nodz(1,1),1); if isempty(column1), column1 = auxMat(auxMat(:,3) == nodz(1,1),1); end; 
    if isempty(column1) == 0, columnWidth1 = columnsDC3(columnsDC3(:,1) == column1(1,1),2); end;
    column2 = auxMat(auxMat(:,2) == nodz(1,2),1); if isempty(column2), column2 = auxMat(auxMat(:,3) == nodz(1,2),1); end; 
    if isempty(column2) == 0, columnWidth2 = columnsDC3(columnsDC3(:,1) == column2(1,1),2); end;
    if (isempty(columnWidth1) & isempty(columnWidth2))
        interiorBeams = [interiorBeams; noFlangeBeams(i,1)];
        continue;
    end
    minColWid = min([columnWidth1, columnWidth2]);
    
    reinfArea = max(beamsDC3(beamsDC3(:,1) == noFlangeBeams(i,1),6));
    b = max(beamsDC3(beamsDC3(:,1) == noFlangeBeams(i,1),3)); %%%%%%%%%%%
    h = max(beamsDC3(beamsDC3(:,1) == noFlangeBeams(i,1),2));%%%%%%%%%%%%%
    d = h - (cover + .02);
    
    width = minColWid;
    
    flangeW = width ; %- max(beamsDC3(beamsDC3(:,1) == beamsBef(i,1),3)); %%%%%%%%%%%%%%
    ratio = reinfArea / (reinfArea + hfAs * flangeW) ;
    reinfPerc = reinfArea * fyd / ( b * d * fcd);
    redBenMom = abaco(ratio, reinfPerc);
    MRbi = redBenMom * b * d^2 * fcd * 1000;
    %%%%%%%%%%%%%%%%%%%%%%%
    direc = element(element(:,1) == noFlangeBeams(i,1),4);
    nodz = fakeBeams(fakeBeams(:,1) == noFlangeBeams(i,1),[2 3]);
    
    auxVec = [];
    for k = 1 : 2
        SmRb = sum(bendingNodes(bendingNodes(:,1) == nodz(1,k) & bendingNodes(:,1) == direc,3));
        SmRc = sum(bendingNodes(bendingNodes(:,1) == nodz(1,k) & bendingNodes(:,1) == 3,3));
        if SmRb > SmRc
            M_Ed = gRd * MRbi * SmRc / SmRc;
        else
            M_Ed = gRd * MRbi;
        end
        auxVec = [auxVec, M_Ed];
    end
    length = fakeBeams(fakeBeams(:,1) == noFlangeBeams(i,1),4);
    z1 = (auxVec(1) - auxVec(2) + DataDesign(DataDesign(:,1) == noFlangeBeams(i,1),4)) / length;
    z2 = (auxVec(2) - auxVec(1) + DataDesign(DataDesign(:,1) == noFlangeBeams(i,1),4)) / length;
    shearForce = max([z1 z2]);
    fakeBeams(fakeBeams(:,1) == noFlangeBeams(i,1),5) = shearForce ;
    
end

%%
for i = 1 : size(fakeBeams,1)
    trueFz = max(fakeBeams(i, 5), DataDesign(DataDesign(:,1) == fakeBeams(i, 1),4));
    b_input = max(beamsDC3(beamsDC3(:,1) == fakeBeams(i, 1),3));
    h_input = max(beamsDC3(beamsDC3(:,1) == fakeBeams(i, 1),2));
    longPhi = max(beamsDC3(beamsDC3(:,1) == fakeBeams(i, 1),5));
    [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = beamDesignDC3shear(fck, fyk , cover, trueFz, b_input, h_input, longPhi);
    beamsDC3(beamsDC3(:,1) == fakeBeams(i,1), [8:11]) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
    
end
%%

for i = 1 : size(barsOfColumns,1)
    for j = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,j) ~= 0
            Fz_Ed = max (DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),3), DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),4));
            b_input = columnsDC3(columnsDC3(:,1) == barsOfColumns(i,jj),2);
            h_input = columnsDC3(columnsDC3(:,1) == barsOfColumns(i,jj),3);
            longPhi = columnsDC3(columnsDC3(:,1) == barsOfColumns(i,jj),5);
            [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = columnDesignDC3shear(fck, fyk , cover, Fz_Ed, b_input, h_input, longPhi);
            
            columnsDC3(columnsDC3(:,1) == barsOfColumns(i,j), [7:10]) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea];
        end
    end
end

%%
origSize = size(beamsDC3,1);
beamsDC3 = [beamsDC3; [barsOfBeams(:,1), zeros(size(barsOfBeams,1) , size(beamsDC3,2)-1)]] ;
times = 1 ;
while times <= noTimesNaming
    for i = (origSize + 1) : size(beamsDC3,1)
        parentBar = barsOfBeams(barsOfBeams(:,1) == beamsDC3(i,1) , 2) ;
        index = find(beamsDC3(:,1) == parentBar) ;
        beamsDC3(i, [2: size(beamsDC3,2)]) = beamsDC3(index, [2: size(beamsDC3,2)]);
        
    end
    times= times + 1;
end

%%%%%%%%%%%clear
toc




%%

auxColumns = unique (columns(:,[2:size(columns,2)]), 'rows');
sizeColumns = size(columns,2);
for i = 1 : size(columns,1)
    for j = 1 : size(auxColumns,1)
        if columns(i, [2: sizeColumns]) == auxColumns(j, [1:size(auxColumns,2)])
            columns(i, sizeColumns + 1) = j ;
        end
    end
end

for p = 1 : size(columns,1)
    columns(p, [8 9]) = element(columns(p,1) == element(:,1),[2 3])
end
