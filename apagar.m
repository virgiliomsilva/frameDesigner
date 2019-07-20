fck = 30;
fyk = 400;
cover = .035;
sec_b = .3;
areaRebar = .5;
longReinfN = 12;
longReinfPh = 20;
[~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, Vrd, sCondition] = DC3columnDesign(fck, fyk , cover, 10, 10, 10, sec_b, [longReinfPh,longReinfN,areaRebar], [], sec_b, longReinfN, longReinfPh)