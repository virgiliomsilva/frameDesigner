% abacusC12_50S500thirdSixth = importdata('info\abacusC12_50S500thirdSixth.csv');
% abacusC12_50S500thirdSixth = scatteredInterpolant(abacusC12_50S500thirdSixth(:,[1:3]),abacusC12_50S500thirdSixth(:,4));
% save('info\abacusC12_50S500thirdSixth.mat','abacusC12_50S500thirdSixth');
% v = .4
% data =  importdata('info\abacusC12_50S500thirdSixth.csv');
% scatter3(data(data(:,1) == v,2), data(data(:,1) == v,3), data(data(:,1) == v,4))

% AA | u | w
%%
% abacusC12_50S500all = importdata('info\abacusC12_50S500all.csv');
% abacusC12_50S500all = scatteredInterpolant(abacusC12_50S500all(:,[1 3]),abacusC12_50S500all(:,2));
% save('info\abacusC12_50S500all.mat','abacusC12_50S500all');

% v = .4
% data =  importdata('info\abacusC12_50S500all.csv');
% scatter3(data(:,1), data(:,3), data(:,2))

%%
% m2 = importdata('beamsEC2_2.mat') ;
% m1 = importdata('beamsEC2_1.mat') ;
% HH = [m1; m2];
% 
% 
% auxxx = [];
% mati = unique (HH(:,1));
% mati(:,[2:11]) = 0;
% for t = 1 : size(mati,1)% unique(HH(:,1))
%     auxxx = HH(HH(:,1) == HH(t,1),:);
%     [val ind] = max(auxxx(:,7));
%     mati(t,:) = auxxx(ind,:);
% end
% mati = sortrows(mati);

%%

% auxColumns = unique (columnsEC2(:,[2:size(columnsEC2,2)]), 'rows');
% sizeColumns = size(columnsEC2,2);
% for i = 1 : size(columnsEC2,1)
%     for j = 1 : size(auxColumns,1)
%         if columnsEC2(i, [2: sizeColumns]) == auxColumns(j, [1:size(auxColumns,2)])
%             columnsEC2(i, sizeColumns + 1) = j ;
%         end
%     end
% end
% 
% for p = 1 : size(columnsEC2,1)
%     columnsEC2(p, [10 11]) = element(columnsEC2(p,1) == element(:,1),[2 3])
% end

%%
axColumns = unique(columnsEC2(:, [2:11]),'rows');
for k = 1 : size(axColumns,1)
    axColumns(k,:) =
    
    
for i = 1 : size(columnsEC2,1)
    for j = 1 : size(axColumns,1)
        