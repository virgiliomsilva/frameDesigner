function [BAR, CASE, FX, FY, FZ, MY, MZ, LENGTH] = infoExtractor(mixed)
trimmed = erase(mixed, '"');
x = split(trimmed, ',');

BarNodeCase = sscanf(string(x(1,1)), '%3d/ %3d/ %3d');

BAR = BarNodeCase(1,1);
CASE = BarNodeCase(3,1);
FX = str2num(cell2mat(x(2,1)));
FY = str2num(cell2mat(x(3,1)));
FZ = str2num(cell2mat(x(4,1)));
MY = str2num(cell2mat(x(5,1)));
MZ = str2num(cell2mat(x(6,1)));
LENGTH = str2num(cell2mat(x(7,1)));