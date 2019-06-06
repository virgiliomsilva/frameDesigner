% example
clear
clc
tic
fnData = 'data\calcada_da_tapada\datasetall.csv' ;
fnNodes = 'data\calcada_da_tapada\nodes.csv' ;
fnElement = 'data\calcada_da_tapada\connectivity.csv' ;

fck = 25 ;
fyk = 500 ;
cover = .035 ;
fyd = fyk / 1.15;
fcd = fck / 1.5;
fywd = fyd;

[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, DataDesign, element, noTimesNaming, stories, nodes, cases] = dataTransformerNew (fnData, fnElement, fnNodes);

beamsEC2 = [];

for bar = fakeBeams([1:50],1)'
    for i = 1 : size(DataDesign,3)

        %bar = 108;
        M_Ed = DataDesign(DataDesign(:,1) == bar,5,i) ;
        Fz_Ed = DataDesign(DataDesign(:,1) == bar,4,i) ;
        length = fakeBeams(fakeBeams(:,1) == bar,4);
        [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = beamDesignEC2new(fck, fyk , cover, M_Ed, Fz_Ed);
        beamNo = bar ; %beamsOnBeams(i,j);
        beamsEC2 = [beamsEC2; [beamNo, sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd,shearReinfPhi, shearReinfSpac, shearReinfLoops]];

    end
end

toc