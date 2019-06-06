function [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, M_Rd] = beamDesignDC3long(fck, fyk , cover, M_Ed, Fz_Ed, length, h_min)
abaco = importdata('info\abacusC12_50S500A1.csv');
longReinforce = importdata('info\steel_beam.csv');

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
%% 

if nargin == 6
    h = floor(length / 12 * 20 ) / 20;
    b = floor(h / 2 * 20) / 20;
elseif nargin == 7
    h = max([floor(length / 12 * 20 ) / 20, h_min]);
    b = floor(h / 2 * 20) / 20;
end

start_h = h;
start_b = b;

roMin = .5 * fctm / fyk;

ratio = 2;
while ratio > 1.05 %%%%%
    %long rebar starting values
    d = h - (cover + .02);  %approximated
    redBendMom = M_Ed / (b * d^2 * fcd * 1000);
    reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
    reinfAreaAux = reinfPerc * b * d * fcd / fyd;
    reinfArea = max([reinfAreaAux, .26*fctm*b*d/fyk, .0013*b*d, roMin*b*h]);
    
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
    h = ceil((final_h * .2 + start_h * .8) * 20) / 20 ;
    start_h = h ;

end