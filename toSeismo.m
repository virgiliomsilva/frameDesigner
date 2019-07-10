function [] = toSeismo(columns, beams, nodes, element, stories, fck, fyk, cover, folder)
%% sections init
sections = {};
%% type columns
columnTypesAux = unique (columns(:,[2:8]),'rows');

infrmFBPH_column = {};
infrmFB_column = {};

for i = 1 : size(columnTypesAux,1)
    infrmFBPH_column(i, [1 : 6]) = {"Column" + i, "Column" + i, 150, 16.67, "None", "0.00"};
    infrmFB_column(i, [1 : 6]) = {"Column" + i, "Column" + i, 5, 150, "None", "0.00"};
    
    h = columnTypesAux(i,1);
    b = columnTypesAux(i,2);
    rebarQtd = columnTypesAux(i,3);
    rebarPhi = columnTypesAux(i,4);
    stirPhi = columnTypesAux(i,5);
    stirSpac = columnTypesAux(i,6);
    stirBran = columnTypesAux(i,7);
    sections{end+1} = sectionWriter(["Column" + i], fyk, fck, h, b, cover, rebarPhi, rebarQtd, stirPhi, stirBran, stirSpac);
end

writetable(cell2table(infrmFBPH_column), [folder '\02_FBPH_columns.csv'],'WriteVariableNames',false);
writetable(cell2table(infrmFB_column), [folder '\02_FB_columns.csv'],'WriteVariableNames',false);
%% type beams
beamTypesAux = unique (beams(:,[2:8]),'rows');

infrmFBPH_beam = {};
infrmFB_beam = {};

for i = 1 : size(beamTypesAux,1)
    infrmFBPH_beam(i, [1 : 6]) = {"Beam" + i, "Beam" + i, 150, 16.67, "None", "0.00"};
    infrmFB_beam(i, [1 : 6]) = {"Beam" + i, "Beam" + i, 5, 150, "None", "0.00"};
    
    h = beamTypesAux(i,1);
    b = beamTypesAux(i,2);
    rebarQtd = beamTypesAux(i,3);
    rebarPhi = beamTypesAux(i,4);
    stirPhi = beamTypesAux(i,5);
    stirSpac = beamTypesAux(i,6);
    stirBran = beamTypesAux(i,7);
    sections{end+1} = sectionWriter(["Beam" + i], fyk, fck, h, b, cover, rebarPhi, rebarQtd, stirPhi, stirBran, stirSpac);
end

writetable(cell2table(infrmFBPH_beam), [folder '\02_FBPH_beams.csv'],'WriteVariableNames',false);
writetable(cell2table(infrmFB_beam), [folder '\02_FB_beams.csv'],'WriteVariableNames',false);
%% sections fin
writetable(cell2table(sections'), [folder '\01_sections.csv'],'WriteVariableNames',false);
%% nodes
for i = 1 : size(nodes,1)
    nodesSeismo(i,[1:5]) = {nodes(i, 1), nodes(i, 2), nodes(i, 3), nodes(i, 4), "structural"};
end

writetable(cell2table(nodesSeismo), [folder '\03_nodes.csv'],'WriteVariableNames',false);
%% element connectivity
for i = 1 : size(element,1)
    if element(i,4) == 3
        [q, sec] = ismember(columns(columns(:,1) == element(i,1), [2:8]), columnTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Column" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    else
        [q, sec] = ismember(beams(beams(:,1) == element(i,1), [2:8]), beamTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Beam" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    end
end

writetable(cell2table(elemConnectivity), [folder '\04_elemConnectivity.csv'],'WriteVariableNames',false);
%% constrains
constrains = {};
for i = setdiff(stories(:,1),0)'
    matAux = nodes(nodes(:,5) == i, 1);
    masterNode = matAux(1);
    slaveNodes = "" ;
    for j = 2 : size(matAux,1)
        slaveNodes = strcat(slaveNodes, "   ",num2str(matAux(j)));
    end
    constrains(end+1, [1:4]) = {"Rigid Diaphragm", masterNode, "X-Y plane", slaveNodes} ;
end

writetable(cell2table(constrains), [folder '\05_constrains.csv'],'WriteVariableNames',false);
%% restraints
restraints = {};
for i = 1 : size(nodes,1)
    if nodes(i,5) == 0
        restraints(end+1,[1 2]) = {nodes(i,1), 'x+y+z+rx+ry+rz'};
    end
end

writetable(cell2table(restraints), [folder '\06_restraints.csv'],'WriteVariableNames',false);