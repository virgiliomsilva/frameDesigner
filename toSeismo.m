%% TO SEISMO
tic
columns = importdata('columnsEC2.mat');
beams = importdata('beamsEC2.mat');
%nodes = importdata('data\nodes.csv');
%element = importdata('data\connectivity.csv');
%% type columns
columnTypesAux = unique (columns(:,[2:9]),'rows');
columnTypes = num2cell(columnTypesAux);

for i = 1 : size(columnTypes,1)
    columnTypes(i, [10 : 15]) = {"Col" + i, "Col" + i, 5, 150, "None", "0.00"};
    columnTypes(i, [16 : 21]) = {"Col" + i, "Col" + i, 150, 16.67, "None", "0.00"};
end

writetable(cell2table(columnTypes), 'toSeismo\columnTypes.csv');

%% type beams
beamTypesAux = unique (beams(:,[2:11]),'rows');
beamTypes = num2cell(beamTypesAux);

for i = 1 : size(beamTypes,1)
    beamTypes(i, [12 : 17]) = {"Beam" + i, "Beam" + i, 5, 150, "None", "0.00"};
    beamTypes(i, [18 : 23]) = {"Beam" + i, "Beam" + i, 150, 16.67, "None", "0.00"};
end

writetable(cell2table(beamTypes), 'toSeismo\beamTypes.csv');

%% nodes
for i = 1 : size(nodes,1)
    nodesSeismo(i,[1:5]) = {nodes(i, 1), nodes(i, 2), nodes(i, 3), nodes(i, 4), "structural"};
end

writetable(cell2table(nodesSeismo), 'toSeismo\nodesSeismo.csv');

%% element connectivity

for i = 1 : size(element,1)
    if element(i,4) == 3
        [q, sec] = ismember(columns(columns(:,1) == element(i,1), [2:9]), columnTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Col" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    else
        [q, sec] = ismember(beams(beams(:,1) == element(i,1), [2:11]), beamTypesAux, 'rows');
        px = element(i,2);
        py = element(i,3);
        elemConnectivity(i, [1:6]) = {element(i,1), "Beam" + sec, px + "   " + py + "   " + "deg=0.00","0.00   0.00   0.00   0.00   0.00   0.00", [],"-1e20   1e20"};
    end
end

writetable(cell2table(elemConnectivity), 'toSeismo\elemConnectivity.csv');

%% constrains & restraints


toc