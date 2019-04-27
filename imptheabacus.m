abacusC12_50S500evenDistr = importdata('info\abacusC12_50S500evenDistr.csv');
abacusC12_50S500evenDistr = scatteredInterpolant(abacusC12_50S500evenDistr(:,[1:3]),abacusC12_50S500evenDistr(:,4));
save('info\abacusC12_50S500evenDistr.mat','abacusC12_50S500evenDistr');