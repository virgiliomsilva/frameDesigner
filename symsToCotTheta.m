% from two syms to one double between 1 and 2.5
function cotTheta = symsToCotTheta(input)
nums = double(input);
vals = nums(nums > 1 & nums < 2.5);

if length(vals) > 1
    cotTheta = max(vals);
else %isempty
    cotTheta = 1;
end