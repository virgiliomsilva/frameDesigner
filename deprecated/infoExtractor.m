function [BAR, CASE, FX, FY, FZ, MY, MZ, LENGTH] = infoExtractor(mixed)
trimmed = erase(mixed, '"');
x = split(trimmed, ',');
%x = erase(x, "'")

BarNodeCase = sscanf(string(x(1,1)), '%3d/ %3d/ %3d');

BAR = BarNodeCase(1,1);
CASE = BarNodeCase(3,1);
FX = str2num(cell2mat(x(2,1)));
FY = abs(str2num(cell2mat(x(3,1))));
FZ = abs(str2num(cell2mat(x(4,1))));
MY = abs(str2num(cell2mat(x(5,1))));
MZ = abs(str2num(cell2mat(x(6,1))));
LENGTH = str2num(cell2mat(x(7,1)));