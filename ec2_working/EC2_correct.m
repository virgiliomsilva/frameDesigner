%% DC2 design for first model
clear
clc

[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories, nodes] = dataTransformer ();

fck = 25 ;
fyk = 500 ;
fyd = fyk / 1.15;
fcd = fck / 1.5;
fywd = fyd;
cover = .035 ;


beamsEC2=[];
for i = size(beamsOnBeams,1) : -1 : 1
    for j = size(beamsOnBeams,2) : -1 : 1
        if beamsOnBeams(i,j) ~= 0
            columnMarker = find(beamsOnBeams(i,:)~=0,1,'last');
            if j == columnMarker
                M_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),5) ;
                Fz_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),4) ;
                length = fakeBeams(fakeBeams(:,1) == beamsOnBeams(i,j),4);
                [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd,shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = beamDesignEC2new(fck, fyk , cover, M_Ed, Fz_Ed);
                savedDim = sec_h;
                beamNo = beamsOnBeams(i,j);
                beamsEC2 = [beamsEC2; [beamNo, sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd,shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd]];
            else
                M_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),5) ;
                Fz_Ed = DataDesign(DataDesign(:,1) == beamsOnBeams(i,j),4) ;
                length = fakeBeams(fakeBeams(:,1) == beamsOnBeams(i,j),4);
                [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd,shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = beamDesignEC2new(fck, fyk , cover, M_Ed, Fz_Ed);
                beamNo = beamsOnBeams(i,j);
                beamsEC2 = [beamsEC2; [beamNo, sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd,shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd]];
            end
        end
    end
    clear columnMarker savedDim 
end

beamsEC2 = unique(beamsEC2,'rows');
%%%%%%%%
origSize = size(beamsEC2,1);
beamsEC2 = [beamsEC2; [barsOfBeams(:,1), zeros(size(barsOfBeams,1) , size(beamsEC2,2)-1)]] ;
times = 1 ;
while times <= noTimesNaming
    for i = (origSize + 1) : size(beamsEC2,1)
        parentBar = barsOfBeams(barsOfBeams(:,1) == beamsEC2(i,1) , 2) ;
        index = find(beamsEC2(:,1) == parentBar) ;
        beamsEC2(i, [2: size(beamsEC2,2)]) = beamsEC2(index, [2: size(beamsEC2,2)]);
        
    end
    times= times + 1;
end

%%%%%%%%
beamsEC2 = unique(beamsEC2,'rows');
save('beamsEC2.mat','beamsEC2')