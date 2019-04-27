%BEAM DESIGN
clear
%clc
tic
%% import data and inputs

%data = 'dataBeams.csv' ; %input('Filename: ');%%%%%%%%%%%%%%%%%
unsData = importdata('data\dataBeams.csv');
abaco = importdata('info\abacusC12_50S500A1.csv');
longReinforce = importdata('info\steel_beam.csv');
shearReinforce = importdata('info\steel_shear.csv');

fck = 25;  %MPa
fyk = 500;  %MPa
cover = .035; 

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd; 
%% transform data
[noRow, noCol] = size(unsData.textdata);
barNameAux = extractBefore(unsData.textdata(2:noRow,:),"/");

for i = 1 : (noRow - 1)
    auxMat(i, 1) = str2double(barNameAux(i, 1));
end

sortData = abs([auxMat, unsData.data]);

row = 1;
rowz = 1;
for j = unique(sortData(:, 1)).'
    for i = 1 : (noRow - 1) 
        if j == sortData(i, 1)
        auxBlock(row, :) = sortData(i, :);
        row = row + 1;
        end
    end
row = 1;
auxMax = max(auxBlock, [], 1);
finalData(rowz, :) = auxMax;
rowz = rowz + 1;
end

noBars = rowz - 1 ;
reinfSolu = zeros(noBars, 29);
clear rowNum colNum unsData barNameAux auxMat sortData row auxBlock auxMax i j rowz data
%% EC 2 design

for i = 1 : noBars
    reinfSolu(i,1) = finalData(i,1);
        
    %long rebar
    M_Ed = finalData(i,5);
    h = .25;
    start_h = h;
    ratio = 2;
    while ratio > 1.05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = floor(.6 * h * 20) / 20;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = h - (cover + .02);  %approximated
        
        %long rebar starting values
        redBendMom = M_Ed / (b * d^2 * fcd * 1000);
        reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
        reinfAreaAux = reinfPerc * b * d * fcd / fyd;
        reinfArea = max([reinfAreaAux, .26 * fctm * b * d / fyk, .0013 * b * d]);
        
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfArea > 0
                diffAuxL(j,1) = longReinforce(j,3) - reinfArea;
            else
                diffAuxL(j,1) = 1000;
            end
        end
            [minDiff, minIndexA] = min(diffAuxL);
            noRebar = longReinforce(minIndexA,2);
            phiRebar = longReinforce(minIndexA,1);
        
        
        %long rebar iterations
        clearance = min([3 , (noRebar - 1)]) * .05;
        M_Rd = 0;
        if M_Ed > M_Rd
            for j = 1 : size(longReinforce,1) %doesn't make sense on the first while iteration, but it comes to use from second on
                if longReinforce(j,3) - reinfArea > 0
                    diffAuxL(j,1) = longReinforce(j,3) - reinfArea;
                else
                    diffAuxL(j,1) = 1000;
                end
                [minDiff, minIndex] = min(diffAuxL);
                reinfSolu(i,4) = longReinforce(minIndex,2);
                reinfSolu(i,5) = longReinforce(minIndex,1);
                reinfSolu(i,6) = longReinforce(minIndex,3);
                %check asmax
            end
           
            if reinfSolu(i,4) <= 4
                noRebar = reinfSolu(i,4);
                space = noRebar - 1;
                clearance = noRebar * .04;
            else
                noRebar = reinfSolu(i,4);
                space = ceil(noRebar / 2) - 1;
                clearance = space * .04;
            end
            
            diff = b - 2 * cover - clearance - noRebar * (phiRebar / 1000);
            while diff < 0
                 b = b + .05;
                 diff = b - 2 * cover - clearance - .02 - noRebar * (phiRebar / 1000); %* reinfSol(i,8) 
            end
            
            bOh = b / h;
            while bOh > .8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                h = h + .05;
                bOh = b / h;
            end
            
            d = h - (cover + .02);
            reinfPerc = reinfSolu(i, 6) * fyd / (b * d * fcd);
            redBendMom = interp1(abaco(:,2).', abaco(:,1).', max([.004, reinfPerc]));
            M_Rd = redBendMom * b * d^2 * fcd *1000;
            
            redBendMom = M_Ed / (b * d^2 * fcd * 1000);
            reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
            reinfAreaAux = reinfPerc * b * d * fcd / fyd;
            reinfArea = max([reinfAreaAux, .26 * fctm * b * d / fyk, .0013 * b * d]);
            
        end
        
        final_h = h;
        ratio = final_h / start_h;
        reinfSolu(i,2) = h;
        reinfSolu(i,3) = b;
        h = ceil((final_h * .2 + start_h * .8) * 20) / 20;
        start_h = h;
        
    end
    reinfSolu(i, 11) = reinfSolu(i, 6) * (reinfSolu(i, 2) - 2 * (cover + .02)) * fyd * 1000;
    
    %stirrups
    z = h - 2 * (cover + .02);  %approximated 
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5)); %p.100 %9.4 '+' 9.5 EC2 %assuming cot(theta) = 2.5
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 %& (mod(shearReinforce(k,2),2) == mod(space + 1,2) | mod(shearReinforce(k,2) / 2,2) == mod(space + 1,2))
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    [minDiffS, minIndexS] = min(diffAuxS,[],1);
    reinfSolu(i,7) = shearReinforce(minIndexS, 1);
    reinfSolu(i,8) = shearReinforce(minIndexS, 3);
    reinfSolu(i,9) = shearReinforce(minIndexS, 2);
    reinfSolu(i,10) = shearReinforce(minIndexS, 4);
    
    redAxial = finalData(i, 2) / (reinfSolu(i, 2) * reinfSolu(i, 3) * fcd * 1000);
    if redAxial > .1
        reinfSolu(i,:) = 666; % string('Not a beam');
        reinfSolu(i,1) = finalData(i,1);
    end
    
end


clear Asw_s b bOh clearance d diff diffAuxL diffAuxS final_h h i j k M_Ed ...
    M_Rd minDiff minDiffS minIndex minIndexA minIndexS noRebar phiRebar ...
    ratio redAxial redBendMom reinfArea reinfAreaAux reinfPerc rowz space start_h z abaco
% adicionar estribos para garantir amarração de 150mm
%SB7zr4RPp
%% EC8 design

roMin = .5 * fctm / fyk ;
sAvaila = unique(shearReinforce(:,3));
for i = 1: noBars
    %long rebar
    reinfAreaEC8 = roMin * reinfSolu(i,2) * reinfSolu(i,3);
    if reinfAreaEC8 > reinfSolu(i,6)
        for j = 1 : size(longReinforce,1)
            if longReinforce(j,3) - reinfAreaEC8 > 0
                diffAux(j,1) = longReinforce(j,3) - reinfAreaEC8;
            else
                diffAux(j,1) = 1000;
            end
        end
        [minDiff, minIndexA] = min(diffAux);
        
        if longReinforce(minIndexA,3) > reinfSolu(i, 6)
            reinfSolu(i, 12) = longReinforce(minIndexA, 2);
            reinfSolu(i, 13) = longReinforce(minIndexA, 1);
            reinfSolu(i, 14) = longReinforce(minIndexA, 3);
        else
            reinfSolu(i, 12) = reinfSolu(i, 4);
            reinfSolu(i, 13) = reinfSolu(i, 5);
            reinfSolu(i, 14) = reinfSolu(i, 6);
        end
    else
        reinfSolu(i, 12) = reinfSolu(i, 4);
        reinfSolu(i, 13) = reinfSolu(i, 5);
        reinfSolu(i, 14) = reinfSolu(i, 6);
    end
    reinfSolu(i, 15) = reinfSolu(i, 14)/ reinfSolu(i, 6) - 1; %increase of steel area in %
    
    %stirrups DC3
    sMatrix = [0 reinfSolu(i, 8); 1 reinfSolu(i, 2) / 4; 2 reinfSolu(i, 7) * .024; 3 .008 * reinfSolu(i, 5)];
    [sMin sIndex] = min(sMatrix(:,2));
    reinfSolu(i, 21) = sMatrix(sIndex, 1);
    for h = 1 : 7
        if sAvaila(h) - sMatrix(sIndex,2) <= 0
            anotherMatrix(i,1) = sAvaila(h) - sMatrix(sIndex,2);
        else
            anotherMatrix(i,1) = -1000;
        end
    end
    
    [sMaxReal sIndexReal] = max(anotherMatrix);
    
    z = reinfSolu(i, 2) - 2 * (cover + .02);
    sw = sAvaila(sIndexReal);
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5));
        
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 & shearReinforce(k,3) == sw
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    
    [minDiffS, minIndexS] = min(diffAuxS);
    reinfSolu(i,17) = shearReinforce(minIndexS, 1);
    reinfSolu(i,18) = shearReinforce(minIndexS, 3);
    reinfSolu(i,19) = shearReinforce(minIndexS, 2);
    reinfSolu(i,20) = shearReinforce(minIndexS, 4);
    reinfSolu(i,22) = reinfSolu(i,20) / reinfSolu(i,10) - 1;
    
  
    %stirrups DC2
    sMatrix = [0 reinfSolu(i, 8); 1 reinfSolu(i, 2) / 4; 2 reinfSolu(i, 7) * .03; 3 .012 * reinfSolu(i, 5)];
    [sMin sIndex] = min(sMatrix(:,2));
    reinfSolu(i, 28) = sMatrix(sIndex, 1);
    for h = 1 : 7
        if sAvaila(h) - sMatrix(sIndex,2) <= 0
            anotherMatrix(i,1) = sAvaila(h) - sMatrix(sIndex,2);
        else
            anotherMatrix(i,1) = -1000;
        end
    end
    
    [sMaxReal sIndexReal] = max(anotherMatrix);
    
    z = reinfSolu(i, 2) - 2 * (cover + .02);
    sw = sAvaila(sIndexReal);
    Asw_s = max(reinfSolu(i,3) * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5));
        
    for k = 1 : size(shearReinforce,1)
        if shearReinforce(k,4) - Asw_s > 0 & shearReinforce(k,3) == sw
            diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
        else
            diffAuxS(k,1) = 1000;
        end
    end
    
    [minDiffS, minIndexS] = min(diffAuxS);
    reinfSolu(i,24) = shearReinforce(minIndexS, 1);
    reinfSolu(i,25) = shearReinforce(minIndexS, 3);
    reinfSolu(i,26) = shearReinforce(minIndexS, 2);
    reinfSolu(i,27) = shearReinforce(minIndexS, 4);
    reinfSolu(i,29) = reinfSolu(i,27) / reinfSolu(i,10) - 1;
end

clear anotherAux auxIndex auxMin cover diffAuxL fcd fck fctm fyd fyk fywd...
    i j k longReinforce minDiff minIndexA noBars noCol noRow reinfAreaEC8...
    roMin sAvaila sIndex sMatrix sMatrix2 sMin shearReinforce diffAux diffAuxS...
    h anotherMatrix minDiffS minIndexS sIndexReal sMaxReal sw z Asw_s

% 
% toc