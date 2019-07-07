%% MrdColumn
% Given a column section it and a axial force it returns the maximum
% uniaxial bending moment
function [M_Rd] = MrdColumn(fck, fyk, b, h, areaRebar, N_Axial)
    if (fyk == 500 && fck > 12 && fck < 50)
        abaco = importdata('info\abacus_REBAP83_C12_C50_S500.mat');
    elseif (fyk == 400 && fck > 12 && fck < 50)
        abaco = importdata('info\abacus_REBAP83_C12_C50_S400.mat');
    else
        error('Materials pair not supported!')
    end

    fcd = fck / 1.5;
    fyd = fyk / 1.15;

    reinfPercFin = areaRebar * fyd / (b * h * fcd);
    redAxialFin = N_Axial / (b * h * fcd * 1000);
    
    diff = zeros(101,2); index = 1;
    for i = 0 : .005 : .5
        reinfPerc = abaco(0, redAxialFin, i);
        diff(index, [1 2]) = [i , abs(reinfPercFin - reinfPerc)];
        index = index + 1;
    end

    [~, index1] = min(diff(:, 2));
    redBenMom = diff(index1, 1);
    M_Rd = redBenMom * h * b^2 * fcd * 1000 ;
end