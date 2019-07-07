function [phiRebar, noRebar, areaRebar] = columnRebar(width, cover, reinfArea, givenLong)
    if exist('givenLong', 'var') && ~isempty(givenLong)
        longReinforce = givenLong;
    else
        longReinforce  = importdata('info\steel_columnEC8.csv');
        longReinforce = longReinforce(:,[1:3]);
    end

    dMax = .03;
    reinfArea = .97 * reinfArea;

    for j = 1 : size(longReinforce,1)
        if longReinforce(j,3) - reinfArea > 0
            diffAux(j) = longReinforce(j,3) - reinfArea;
        else
            diffAux(j) = 1000;
        end
    end

    [val, minIndex] = min(diffAux);
    noRebar = longReinforce(minIndex,2);
    phiRebar = longReinforce(minIndex,1);
    areaRebar = longReinforce(minIndex,3);

    diffe = width - 2 * (cover + .01) - (phiRebar/1000 * (noRebar/4 + 1)) - (max([.02, phiRebar/1000, dMax+.005]) * (noRebar/4));

    if val > 1 || diffe < 0
        phiRebar = 0;
        noRebar = 0;
        areaRebar = 0;
    end
end