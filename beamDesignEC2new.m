function [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, V_Rd] = beamDesignEC2new(fck, fyk , cover, M_Ed, Fz_Ed)
abaco = importdata('info\abacusC12_50S500A1.csv');
longReinforce = importdata('info\steel_beam.csv');
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd;
%% 
% 
% if nargin == 6
%     h = floor(length / 12 * 20 ) / 20;
%     b = floor(h / 2 * 20) / 20;
% elseif nargin == 7
%     h = max([floor(length / 12 * 20 ) / 20, h_min]);
%     b = floor(h / 2 * 20) / 20;
% end

h = .2;
b = .2;

start_h = h;
start_b = b;



ratio = 2;
while ratio > 1.05 %%%%%
    %long rebar starting values
    d = h - (cover + .02);  %approximated
    redBendMom = M_Ed / (b * d^2 * fcd * 1000);
    reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
    reinfAreaAux = reinfPerc * b * d * fcd / fyd;
    reinfArea = max([reinfAreaAux, .26*fctm*b*d/fyk, .0013*b*d]);
    
    %long rebar iterations
    M_Rd = 0;
    if M_Ed > M_Rd
        for j = 1 : size(longReinforce,1) 
            if longReinforce(j,3) - reinfArea > 0
                diffAuxL(j,1) = longReinforce(j,3) - reinfArea;
            else
                diffAuxL(j,1) = 1000;
            end
            [minDiff, minIndex] = min(diffAuxL);
            longReinfNo = longReinforce(minIndex,2);
            longReinfPhi = longReinforce(minIndex,1);
            longReinfArea = longReinforce(minIndex,3);
            %check asmax
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if longReinfNo <= 4
            noRebar = longReinfNo;
            space = noRebar - 1;
            clearance = noRebar * .05;
        else
            noRebar = longReinfNo;
            space = ceil(noRebar / 2) - 1;
            clearance = space * .05;
        end

        diff = b - 2 * cover - clearance - longReinfNo * (longReinfPhi / 1000);
        
        while diff < 0
             b = b + .05;
             diff = b - 2 * cover - clearance - .02 - longReinfNo * (longReinfPhi / 1000);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bOh = b / h ;
        while bOh > .8
            h = h + .05;
            bOh = b / h;
        end

        d = h - (cover + .02);
        reinfPerc = longReinfArea * fyd / (b * d * fcd);
        redBendMom = interp1(abaco(:,2).', abaco(:,1).', max([.004, reinfPerc]));
        M_Rd = redBendMom * b * d^2 * fcd *1000 ;
    end

    
    final_h = h ;
%     final_b = b ;
    ratio = final_h / start_h;
    sec_h = h;
    sec_b = b ;
    h = ceil((final_h * .05 + start_h * .95) * 20) / 20 ;
    start_h = h ;
    b = floor(.5 * h * 20) / 20;

end



%stirrups
z = h - 2 * (cover + .02);  %approximated 
Asw_s = max(sec_b * .08 * sqrt(fck) / fyk, Fz_Ed / (z * fywd * 1000 * 2.5)); %p.100 %9.4 '+' 9.5 EC2 %assuming cot(theta) = 2.5
for k = 1 : size(shearReinforce,1)
    if shearReinforce(k,4) - Asw_s > 0 %& (mod(shearReinforce(k,2),2) == mod(space + 1,2) | mod(shearReinforce(k,2) / 2,2) == mod(space + 1,2))
        diffAuxS(k,1) = shearReinforce(k,4) - Asw_s;
    else
        diffAuxS(k,1) = 1000;
    end
end
[minDiffS, minIndexS] = min(diffAuxS,[],1);
shearReinfPhi = shearReinforce(minIndexS, 1);
shearReinfSpac = shearReinforce(minIndexS, 3);
shearReinfLoops = shearReinforce(minIndexS, 2);
shearReinfArea = shearReinforce(minIndexS, 4);
V_Rd = shearReinfArea * z * fywd * 2.5 / shearReinfSpac;