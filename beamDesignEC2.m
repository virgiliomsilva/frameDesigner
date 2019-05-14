function [sec_h, sec_b, longReinfNo, longReinfPhi, longReinfArea, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea] = ...
    beamDesignEC2(fck, fyk , cover, M_Ed, Fz_Ed, b_input, h_input)
%%
% clear
% clc
%%
% fck = 25;
% fyk = 500;
% cover = 0.035;
% M_Ed = 150;
% Fz_Ed = 30;
% b_input = .2;
%%
abaco = importdata('info\abacusC12_50S500A1.csv');
longReinforce = importdata('info\steel_beam.csv');
shearReinforce = importdata('info\steel_shear.csv');

fcd = fck / 1.5; 
fctm = .3 * fck^(2/3);
fyd = fyk / 1.15;
fywd = fyd; 
%% % EC 2 design
%reinfSolu(i,1) = finalData(i,1);

%long rebar
%M_Ed %= finalData(i,5);
if nargin == 7
    h = h_input ;
else
    h = min([ .20, b_input]); % at minimun is a square
end

start_h = h;
ratio = 2;
b = b_input;

while ratio > 1.05 %%%%%
    
    %delete rows of longReinforce that don't fit on the given b_input
    deleteRow = [];
    for ii = 1 : size(longReinforce,1) 
        if longReinforce(ii, 2) <= 3
            noBarras = longReinforce(ii, 2);
            espaco = noBarras - 1 ;
            espacoEntre = espaco * .05 ;
        else
            noBarras = longReinforce(ii, 2) ;
            espaco = ceil(noBarras / 2) - 1 ;
            espacoEntre = espaco * .05 ;
        end
        
        dife = b - 2 * cover - espacoEntre - noBarras * (longReinforce(ii, 1) / 1000);
        
        if dife < 0
            deleteRow = [deleteRow, ii] ;
        end
    end
    longReinforce(deleteRow,:) = [] ;   
    
    %long rebar starting values
    d = h - (cover + .02);  %approximated
    redBendMom = M_Ed / (b * d^2 * fcd * 1000);
    reinfPerc = interp1(abaco(:,1).',abaco(:,2).',max([.005, redBendMom]));
    reinfAreaAux = reinfPerc * b * d * fcd / fyd;
    reinfArea = max([reinfAreaAux, .26 * fctm * b * d / fyk, .0013 * b * d]);
    
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
    ratio = final_h / start_h;
    sec_h = h ;
    sec_b = b ;
    h = ceil((final_h * .2 + start_h * .8) * 20) / 20 ;
    start_h = h ;

end


%col11 = col6 * col2 - 2 * (cover + .02)) * fyd * 1000;


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

% redAxial = finalData(i, 2) / (reinfSolu(i, 2) * reinfSolu(i, 3) * fcd * 1000);
% if redAxial > .1
%     reinfSolu(i,:) = 666; % string('Not a beam');
%     reinfSolu(i,1) = finalData(i,1);
% end
    



% clear Asw_s b bOh clearance d diff diffAuxL diffAuxS final_h h i j k M_Ed ...
%     M_Rd minDiff minDiffS minIndex minIndexA minIndexS noRebar phiRebar ...
%     ratio redAxial redBendMom reinfArea reinfAreaAux reinfPerc rowz space start_h z abaco
% % adicionar estribos para garantir amarração de 150mm
% %SB7zr4RPp


