%% designer 
% function [] = frameDesignerDC1 (fnData, fnElement, fnNodes, fck, fyk, cover)
%%DELETE
clear
clc
tic

% % fnData = 'data\calcada_da_tapada\datasetall.csv' ;
% % fnNodes = 'data\calcada_da_tapada\nodes.csv' ;
% % fnElement = 'data\calcada_da_tapada\connectivity.csv' ;

fnData = 'data\regular\datasetall.csv' ;
fnNodes = 'data\regular\nodes.csv' ;
fnElement = 'data\regular\connectivity.csv' ;

% fck = 30 ;
% fyk = 400 ;
% cover = .035 ;

%%
[barsOfBeams, barsOfColumns, beamDesiOrd, beamsOnBeams, fakeBeams, ...
    DataDesign, element, noTimesNaming, stories, nodes, cases] = ...
    dataTransformerNew (fnData, fnElement, fnNodes);


% fyd = fyk / 1.15;
% fcd = fck / 1.5;
% fywd = fyd;

% while
   
% end

disp(['** Finished in ', num2str(round(toc,2)), ' secs **']);