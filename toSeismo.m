%% TO SEISMO
% function [] = toSeismo(columns, beams, nodes, element, stories)

tic
columns = importdata('columnsEC2.mat');
beams = importdata('beamsEC2.mat');
nodes = importdata('data\regular\nodes.csv');
element = importdata('data\regular\connectivity.csv');
%% type columns
columnTypesAux = unique (columns(:,[2:9]),'rows');

infrmFBPH_column = {};
infrmFB_column = {};

for i = 1 : size(columnTypesAux,1)
    infrmFBPH_column(i, [1 : 6]) = {"Column" + i, "Column" + i, 150, 16.67, "None", "0.00"};
    infrmFB_column(i, [1 : 6]) = {"Column" + i, "Column" + i, 5, 150, "None", "0.00"};
end

writetable(cell2table(infrmFBPH_column), 'toSeismo\02_FBPH_columns.csv','WriteVariableNames',false);
writetable(cell2table(infrmFB_column), 'toSeismo\02_FB_columns.csv','WriteVariableNames',false);
%% type beams
beamTypesAux = unique (beams(:,[2:11]),'rows');

infrmFBPH_beam = {};
infrmFB_beam = {};

for i = 1 : size(beamTypesAux,1)
    infrmFBPH_beam(i, [1 : 6]) = {"Beam" + i, "Beam" + i, 150, 16.67, "None", "0.00"};
    infrmFB_beam(i, [1 : 6]) = {"Beam" + i, "Beam" + i, 5, 150, "None", "0.00"};
end

writetable(cell2table(infrmFBPH_beam), 'toSeismo\02_FBPH_beams.csv','WriteVariableNames',false);
writetable(cell2table(infrmFB_beam), 'toSeismo\02_FB_beams.csv','WriteVariableNames',false);
%% nodes
for i = 1 : size(nodes,1)
    nodesSeismo(i,[1:5]) = {nodes(i, 1), nodes(i, 2), nodes(i, 3), nodes(i, 4), "structural"};
end

writetable(cell2table(nodesSeismo), 'toSeismo\03_nodes.csv','WriteVariableNames',false);
%% element connectivity
for i = 1 : size(element,1)
    if element(i,4) == 3
        [q, sec] = ismember(columns(columns(:,1) == element(i,1), [2:9]), columnTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Column" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    else
        [q, sec] = ismember(beams(beams(:,1) == element(i,1), [2:11]), beamTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Beam" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    end
end

writetable(cell2table(elemConnectivity), 'toSeismo\04_elemConnectivity.csv','WriteVariableNames',false);
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

writetable(cell2table(constrains), 'toSeismo\05_constrains.csv','WriteVariableNames',false);
%% restraints
restraints = {};
for i = 1 : size(nodes,1)
    if nodes(i,5) == 0
        restraints(end+1,[1 2]) = {nodes(i,1), 'x+y+z+rx+ry+rz'};
    end
end

writetable(cell2table(restraints), 'toSeismo\06_restraints.csv','WriteVariableNames',false);
toc