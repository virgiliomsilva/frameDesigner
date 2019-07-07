%% COMPATIBILITY FUNCTION
% Given a steel reinforcement pattern the function returns
% the other possible patterns suitable to merge with the given on.
% Pattern may be inserted as total area or number of long reinforcement
% rods and phi;
% Flag inputs may be 'EC2' or 'EC8'
% Input4 may be inserted to eleiminate possibilities with lesser area
function [newComp] = columnComp(input1, input2, input3, input4)
    if nargin == 2
        %input1 area  %input2 flag
        if strcmpi(input2, 'EC8')
            table = importdata('info\steel_columnEC8.csv');
        elseif strcmpi(input2, 'EC2')
            table = importdata('info\steel_column.csv');
        end

        % get reference pattern
        [row,~] = find(table(:,3) == input1);

    elseif nargin == 3
        %input1 no long reinfo  %input2 phi long reinf %input3 flag
        if strcmpi(input3, 'EC8')
            table = importdata('info\steel_columnEC8.csv');
        elseif strcmpi(input3, 'EC2')
            table = importdata('info\steel_column.csv');
        end

        % get rference pattern
        firstIdx = ismember(table(:,2), input1);
        seconIdx = ismember(table(:,1), input2);
        [row,~]  = find(logical(floor(sum([firstIdx, seconIdx],2) / 2)));

    elseif nargin == 4
        % input1 no long reinfo  %input2 phi long reinf %input3 flag %input4 cutarea
        % if input 4 exists eliminate possible area below that one
        if strcmpi(input3, 'EC8')
            table = importdata('info\steel_columnEC8.csv');
        elseif strcmpi(input3, 'EC2')
            table = importdata('info\steel_column.csv');
        end

        % get rference pattern
        firstIdx = ismember(table(:,2), input1);
        seconIdx = ismember(table(:,1), input2);
        [row,~]  = find(logical(floor(sum([firstIdx, seconIdx],2) / 2)));
        % delete rows with a lower area
        toDel = [];
        for i = 1 : size(table,1)
            if table(i,3) < table(row, 3)
                toDel(1,end+1) = i ;
            end
        end
        table(toDel, :) = [];
        % refind reference pattern row
        firstIdx = ismember(table(:,2), input1);
        seconIdx = ismember(table(:,1), input2);
        [row,~]  = find(logical(floor(sum([firstIdx, seconIdx],2) / 2)));

    else
        error('Refactor the inputs!');
    end

    if length(row) > 1
        error('Ambiguous area! Try using 3 inputs');
    end

    %%
    newComp = [];
    %possible diameters & quantities
    phiTable = [12, 16, 20, 25, 32, 32];
    % qtdTable = [4, 8, 12, 16, 20, 24];
    %%
    % divide in regular and non regular pattern
    nonRegularRowStart = find(table(:,4) > 0,1);
    part1 = table([1: nonRegularRowStart-1], [1:3]);
    part2 = table([nonRegularRowStart:end], :);

    if row >= nonRegularRowStart
        %know what are we working with
        mainPhi = table(row, 5);
        originalQtd = table(row, 2);
        %get which combinations are possible
        if mod(originalQtd/4 + 1,2) == 0
            possQtd = [originalQtd, originalQtd+4];%don't want a big change of quantity
        else
            possQtd = [originalQtd-4, originalQtd];
        end

        if mainPhi>12
            index = find(phiTable == mainPhi, 1, 'first');
            possPhi = phiTable([index-1, index]);
        else
            possPhi = 12;
        end
        %remove possibilites that are not compatible
        %intersects the two previous conditions and gets the rows that follow
        firstIdx = ismember(part1(:,1), possPhi);
        seconIdx = ismember(part1(:,2), possQtd);
        newComp = [newComp; part1(logical(floor(sum([firstIdx, seconIdx],2) / 2)), [1:3])];
        firstIdx = ismember(part2(:,5), mainPhi);
        seconIdx = ismember(part2(:,2), possQtd);%%%%%%modified here
        newComp = [newComp; part2(logical(floor(sum([firstIdx, seconIdx],2) / 2)), [1:3])];
    else
        %know what are we working with
        mainPhi = table(row, 1);
        originalQtd = table(row, 2);
        %get which combinations are possible
        if mod(originalQtd/4 + 1,2) == 0
            possQtd = [originalQtd, originalQtd+4];%don't want a big change of quantity
        else
            possQtd = [originalQtd-4, originalQtd];
        end

        if mainPhi>12
            index = find(phiTable == mainPhi, 1, 'first');
            possPhi = phiTable([index-1, index, index+1]);
        else
            index = find(phiTable == mainPhi);
            possPhi = phiTable([index, index+1]);
        end
        %remove possibilites that are not compatible
        firstIdx = ismember(part1(:,1), possPhi);
        seconIdx = ismember(part1(:,2), possQtd);
        newComp = [newComp; part1(logical(floor(sum([firstIdx, seconIdx],2) / 2)), [1:3])];
        firstIdx = ismember(part2(:,5), phiTable([index, index +1]));
        seconIdx = ismember(part2(:,2), possQtd);%%% also modified here
        newComp = [newComp; part2(logical(floor(sum([firstIdx, seconIdx],2) / 2)), [1:3])];
    end
end