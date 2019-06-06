% % % % info = importdata('Book.csv');
% % % % x1 = info(:,1);
% % % % y1 = info(:,2);
% % % % z1 = info(:,3);
% % % % x2 = info([1:18],6);
% % % % y2 = info([1:18],7);
% % % % z2 = info([1:18],8);
% % % % hold on
% % % % scatter3(x1,y1,z1,'r');
% % % % scatter3(x2,y2,z2,'b');
% % % % hold off
% % % 
% % % % range = [1:11];
% % % % x1 = abacusC12_50S500evenDistr.Points(range,2);
% % % % y1 = abacusC12_50S500evenDistr.Points(range,3);
% % % % z1 = abacusC12_50S500evenDistr.Values(range,1);
% % % % scatter3(x1,y1,z1,'r');
% % % 
% % % allVals = abacusC12_50S500evenDistr.Points;
% % % allVals(:,4) = abacusC12_50S500evenDistr.Values(:,1);
% % % 
% % % w1 = allVals;
% % % w1(w1(:,4) ~= 1 , :) = [];
% % 
% % z1 = w1(:,1);
% % y1 = w1(:,2);
% % x1 = w1(:,3);
% % scatter3(x1,y1,z1,'r');
% 
% % wut = menu('which op?','op1','op2')
% a = 3;
% 
% switch a
%     case 2
%         disp('sda')
%     otherwise disp('other')
% end
  


% tf = strncmpi(s1,s2,4)
sections = {};

for i = 1 : 12: 100
    sections{end+1} = i+1
end