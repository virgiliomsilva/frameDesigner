function [] = beamDesignEC2(M_Ed,)


%%
abaco = importdata('abacusC12_50S500A1.csv');
longReinforce = importdata('steel_beam.csv');
shearReinforce = importdata('steel_shear.csv');

fck = 25;  %MPa
fyk = 500;  %MPa
cover = .035; 

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd; 
%%
% EC 2 design
%reinfSolu(i,1) = finalData(i,1);

%long rebar
M_Ed %= finalData(i,5);
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
            clearance = noRebar * .05;
        else
            noRebar = reinfSolu(i,4);
            space = ceil(noRebar / 2) - 1;
            clearance = space * .05;
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

    end

    final_h = h;
    ratio = final_h / start_h;
    col2 = h;
    col3 = b;
    h = ceil((final_h * .2 + start_h * .8) * 20) / 20;
    start_h = h;

end
col11 = col6 * col2 - 2 * (cover + .02)) * fyd * 1000;

%stirrups
z = h - 2 * (cover + .02);  %approximated 
Asw_s = max(col3 * .08 * sqrt(fck) / fyk, finalData(i, 4) / (z * fywd * 1000 * 2.5)); %p.100 %9.4 '+' 9.5 EC2 %assuming cot(theta) = 2.5
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
    



clear Asw_s b bOh clearance d diff diffAuxL diffAuxS final_h h i j k M_Ed ...
    M_Rd minDiff minDiffS minIndex minIndexA minIndexS noRebar phiRebar ...
    ratio redAxial redBendMom reinfArea reinfAreaAux reinfPerc rowz space start_h z abaco
% adicionar estribos para garantir amarração de 150mm
%SB7zr4RPp
