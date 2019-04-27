%% function caller
clear
clc
tic

[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories] = dataTransformer ();

fck = 25 ;
fyk = 500 ;
cover = .035 ;
%% columns
columns = [];
columnsAux = [];
columnsAux2 = [];
for i = 1 : size(barsOfColumns,1)
    for j = 1 : size(barsOfColumns,2)
        b = columnDesignEC2(fck, fyk, cover, DataDesign(DataDesign(:,1) == barsOfColumns(i,j),2), ...
            DataDesign(DataDesign(:,1) == barsOfColumns(i,j),5), DataDesign(DataDesign(:,1) == barsOfColumns(i,j),6));
        bAux (j, 1) = b;
    end
    bColumn = max(bAux) ;
    
    for jj = 1 : size(barsOfColumns,2)
        bar = barsOfColumns(i, jj);
        [b, noRebar, phiRebar, areaRebar, phiShear, spacingShear] = columnDesignEC2(fck, ...
            fyk, cover, DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),2), ...
            DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),5), ...
            DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),6), bColumn);
        columnsAux = [columnsAux; [bar, b, noRebar, phiRebar, areaRebar, phiShear, spacingShear]] ;
    end
    
    [maxV indexMax] = max(columnsAux(:, 5));
    [minV indexMin] = min(columnsAux(:, 7));
    for k = 1 : size(columnsAux,1)
        columnsAux2(k, [1 2]) = columnsAux(k, [1 2]); 
        columnsAux2(k, [3 4]) = columnsAux(indexMax, [3 4]); 
        columnsAux2(k, [5 6]) = columnsAux(indexMin, [6 7]); 
    end
    
    columns = [columns; columnsAux2];
    clear bColumn bAux %columnsAux columnsAux2
    columnsAux = [];
    columnsAux2 = [];
    
end
columns(columns(:,1) == 0,:) = [] ;
%columns = unique(columns, 'rows');

clear areaRebar b bar columnsAux columnsAux2 i indexMax indexMin j jj k ...
    maxV minV noRebar phiRebar phiShear spacingShear 
%% beams

beams = [];
auxVec1 = [];
auxVec2 = [];
auxVec3 = [];
auxVec4 = [];
auxBlock = element(element(:,4) == 3, :);
for i = 1 : size(beamDesiOrd,2)
    auxVec1 = [auxVec1, beamDesiOrd(1,i)];
    
    for j = 1 : noTimesNaming %bars of i beam
        if j <= size(auxVec1,2)
        iniBar = barsOfBeams(barsOfBeams(:,2) == auxVec1(1,j),1);
        auxVec1 = [auxVec1, iniBar];
        end
    end
    
    for k = 1 : size(auxVec1,2) % nodes of those bars
       nodz = element(element(:,1) == auxVec1(1,k),[2 3]);
       auxVec2 = [auxVec2, nodz] ;
    end
    auxVec2 = auxVec2' ;
    auxVec2 = unique(auxVec2, 'rows') ;
    auxVec2 = auxVec2' ;
    
    for jj = 1 : size(auxVec2,2) % columns with those nodes
        %if isempty(auxBlock(:,[2 3]) == auxVec2(1,jj)) == 1
        iniElem = auxBlock(auxBlock(:,2) == auxVec2(1,jj),1) ;
        auxVec3 = [auxVec3; iniElem] ;
        iniElem = auxBlock(auxBlock(:,3) == auxVec2(1,jj),1) ;
        auxVec3 = [auxVec3;, iniElem] ;
        %end
    end
    auxVec3 = unique(auxVec3, 'rows') ;
    auxVec3 = auxVec3' ;
    
    if isempty(auxVec3) == 1
        min_b = .2 ;
    else    
        for ii = 1 : size(auxVec3,1) % b width of those columns
            AuxB = columns(columns(:,1) == auxVec3(ii,1),2) ;
            auxVec4 = [auxVec4, AuxB] ;
        end

        min_b = min(auxVec4) ;
    end
    auxVec1 = [];
    auxVec2 = [];
    auxVec3 = [];
    auxVec4 = [];
    
    % beam design per se
    bar = beamDesiOrd(1, i) ;
    M_Ed = DataDesign(DataDesign(:,1) == beamDesiOrd(1,i),5) ;
    Fz_Ed = DataDesign(DataDesign(:,1) == beamDesiOrd(1,i),4) ;
    b_input = min_b ;
    
    [col2, col3, col4, col5, col6, col7, col8, col9, col10] = beamDesignEC2(fck, fyk , cover, M_Ed, Fz_Ed, b_input) ;
    
    %beams = [beams; [bar, col2, col3, col4, col5, col6, col7, col8, col9, col10]] ;
    beams = [beams; [bar, col2, col3, col4, col5, col7, col8, col9]] ;
end

origSize = size(beams,1);
beams = [beams; [barsOfBeams(:,1), zeros(size(barsOfBeams,1) , size(beams,2)-1)]] ;
times = 1 ;
while times <= noTimesNaming
    for i = (origSize + 1) : size(beams,1)
        parentBar = barsOfBeams(barsOfBeams(:,1) == beams(i,1) , 2) ;
        index = find(beams(:,1) == parentBar) ;
        beams(i, [2: size(beams,2)]) = beams(index, [2: size(beams,2)]);
        
    end
    times= times + 1;
end

%beams = unique(beams, 'rows');

clear AuxB auxBlock auxVec auxVec2 auxVec3 auxVec4 b_input bar barsOfBeams ...
    barsOfColumns beamsOnBeams col10 col2 col3 col4 col5 col6 ...
    col7 col8 col9 cover fck fyk Fz_Ed i ii index iniBar ...
    iniElem j jj k M_Ed min_b nodz noTimesNaming origSize parentBar stories times ...
    % element DataDesign beamDesiOrd
%%
% auxColumns = unique (columns(:,[2:size(columns,2)]), 'rows');
% sizeColumns = size(columns,2);
% for i = 1 : size(columns,1)
%     for j = 1 : size(auxColumns,1)
%         if columns(i, [2: sizeColumns]) == auxColumns(j, [1:size(auxColumns,2)])
%             columns(i, sizeColumns + 1) = j ;
%         end
%     end
% end
% 
% for p = 1 : size(columns,1)
%     columns(p, [8 9]) = element(columns(p,1) == element(:,1),[2 3])
% end


%%
% auxBeams = unique (beams(:,[2:size(beams,2)]), 'rows');
% sizeBeams = size(beams,2);
% for i = 1 : size(beams,1)
%     for j = 1 : size(auxBeams,1)
%         if beams(i, [2: sizeBeams]) == auxBeams(j, [1:size(auxBeams,2)])
%             beams(i, sizeBeams + 1) = j ;
%         end
%     end
% end
% 
% for p = 1 : size(beams,1)
%     beams(p, [10 11]) = element(beams(p,1) == element(:,1),[2 3])
% end


toc