clear
clc

[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories, nodes] = dataTransformer ();


fck = 25 ;
fyk = 500 ;
fyd = fyk / 1.15;
fcd = fck / 1.5;
fywd = fyd;
cover = .035 ;


barsofcolumns2 = [barsOfColumns(:,1);barsOfColumns(:,2);barsOfColumns(:,3);barsOfColumns(:,4);barsOfColumns(:,5)];
barsofcolumns2(barsofcolumns2(:,1) == 0,:) = [];

beamsEC2 = importdata('beamsEC2.mat');
auxVec1 = [];
auxVec2 = [];
auxVec3 = [];
auxVec4 = [];
auxBlock = element(element(:,4) ~= 3, :);

columnsEC2 = [];
bAux = [];
for i = 1 : size(barsOfColumns,1)
    for j = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,j) ~= 0
%             N_Axial = DataDesign(DataDesign(:,1) == barsOfColumns(i,j),2);
%             My_h = max([DataDesign(DataDesign(:,1) == barsOfColumns(i,j),5)]);
%             Mz_b = max([DataDesign(DataDesign(:,1) == barsOfColumns(i,j),6)]);
%             [h, b, noRebar, phiRebar, areaRebar, M_Rd, phiShear, spacingShear] = columnDesignEC2new(fck, fyk , cover, N_Axial, My_h, Mz_b);
%             bAux = [bAux; [h, b, noRebar, phiRebar]];
            %%%%%%%
            auxVec1 = [auxVec1, barsOfColumns(i,j)];

%             for j = 1 : noTimesNaming %bars of i beam
%                 if j <= size(auxVec1,2)
%                 iniBar = barsOfBeams(barsOfBeams(:,2) == auxVec1(1,j),1);
%                 auxVec1 = [auxVec1, iniBar];
%                 end
%             end

            for k = 1 : size(auxVec1,2) % nodes of those bars
               nodz = element(element(:,1) == auxVec1(1,k),[2 3]);
               auxVec2 = [auxVec2, nodz] ;
            end
            auxVec2 = auxVec2' ;
            auxVec2 = unique(auxVec2, 'rows') ;
            auxVec2 = auxVec2' ;

            for jj = 1 : size(auxVec2,2) % beams with those nodes
                %if isempty(auxBlock(:,[2 3]) == auxVec2(1,jj)) == 1
                iniElem = auxBlock(auxBlock(:,2) == auxVec2(1,jj),1) ;
                auxVec3 = [auxVec3; iniElem] ;
                iniElem = auxBlock(auxBlock(:,3) == auxVec2(1,jj),1) ;
                auxVec3 = [auxVec3; iniElem] ;
                %end
            end
            auxVec3 = unique(auxVec3, 'rows') ;
            auxVec3 = auxVec3' ;

            if isempty(auxVec3)
                max_b = .2 ;
            else    
                for ii = 1 : size(auxVec3,1) % b width of those columns
                    AuxB = beamsEC2(beamsEC2(:,1) == auxVec3(ii,1),3) ;
                    auxVec4 = [auxVec4; AuxB] ;
                end

               % max_b = max(auxVec4) ;
            end
            auxVec1 = [];
            auxVec2 = [];
            auxVec3 = [];
            %auxVec4 = [];
            
            
            %%%%%%%
        end
    end
    
    max_b = max(auxVec4);
    auxVec4 = [];
%     hColumn = max_b;
%     bColumn = max_b;
%     phiColumn = min(bAux(:,4));
    
    for jj = 1 : size(barsOfColumns,2)
        if barsOfColumns(i,jj) ~= 0
            bar = barsOfColumns(i, jj);
            N_Axial = DataDesign(DataDesign(:,1) == barsOfColumns(i,jj),2);
            My_h = max([DataDesign(DataDesign(:,1) == barsOfColumns(i,j),5)]);
            Mz_b = max([DataDesign(DataDesign(:,1) == barsOfColumns(i,j),6)]);
            [h, b, noRebar, phiRebar, areaRebar, M_Rd, phiShear, spacingShear] = columnDesignEC2new(fck, fyk , cover, N_Axial, My_h, Mz_b, max_b);
            columnsEC2 = [columnsEC2; [bar, h, b, noRebar, phiRebar, areaRebar, M_Rd, phiShear, spacingShear]];
            %bendingNodes = [bendingNodes; [element(element(:,1) == barsOfColumns(i,jj),2), 3, M_Rd]; [element(element(:,1) == barsOfColumns(i,jj),3), 3, M_Rd]];
        end
    end
    clear bColumn hColumn phiColumn
    bAux = [];
end

save('columnsEC2.mat','columnsEC2')