auxBlock = element(element(:,4) == 3,:)

%%
auxBlock = [auxBlock, zeros(size(auxBlock,1),4)]

%%
for i = 1 : size(auxBlock,1)
   auxBlock(i,5) = nodes(nodes(:,1) == auxBlock(i,2), 4)
   auxBlock(i,6) = nodes(nodes(:,1) == auxBlock(i,3), 4) 
end

for i = 1 : size(auxBlock,1)
    auxBlock(i,7) = max([auxBlock(i,5) auxBlock(i,6)])
end

auxBlock(auxBlock(:,7) ~= 3.05,:) = [] % piso

for i = 1 : size(auxBlock,1)
    auxBlock(i,8) = DataDesign(DataDesign(:,1) == auxBlock(i,1), 2)
end