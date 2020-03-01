%% DC1frameDesigner
function [] = DC1frameDesigner(buildingName, seismicVerticalLoadCase, nonSeismicCases, folder)
tic; disp('Started DC1');
%% STEP 1 - GETTING DATA
loading = waitbar(0,'Reading data','Name', 'DC1: Step 1 of 4');

fnData    = ['data\' buildingName '\dataset.csv'] ;
fnNodes   = ['data\' buildingName '\nodes.csv'] ;
fnElement = ['data\' buildingName '\connectivity.csv'] ;

[~, barsOfColumns, beamDesiOrd, ~, ~, DataDesignMax, DataDesignMin, element, ~, stories, nodes, cases] = dataTransformer (fnData, fnElement, fnNodes);

clear buildingName fnData fnElement fnNodes

allCasesIdx                = [1:length(cases)];
seismicVerticalLoadCaseIdx = find(ismember(cases, seismicVerticalLoadCase));
nonSeismicCasesIdx         = find(ismember(cases, nonSeismicCases));

loading = waitbar(1,loading,'Reading data','Name', 'DC1: Step 1 of 4'); pause(1);
%% STEP 2 - DESIGN BEAMS
close(loading); loading = waitbar(0,'Initializing beams','Name', 'DC1: Step 2 of 4'); pause(1);

beams = []; beamsMid = [];
for i = 1 : length(beamDesiOrd)
    barIndex = find(DataDesignMax(:,1,1) == beamDesiOrd(i));
    
    %longitudinal rebar and stirrups on critical region
    maxM_Ed = max(abs(DataDesignMax(barIndex, 5, :)));
    Fz_Ed = max(max(abs(DataDesignMax(barIndex, 4, allCasesIdx)),abs(DataDesignMin(barIndex, 4, allCasesIdx))));
    
    [sec_h, sec_b, longReinfNo, longReinfPhi, ~, M_Rd, roMinCondition, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition] = ...
        DC1beamDesign(maxM_Ed, Fz_Ed);
    
    %midStirrups
    [shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid] = DC1beamDesignMidShear(Fz_Ed, sec_b, sec_h, longReinfNo, longReinfPhi);
    
    beams   (end+1,:) = [DataDesignMax(barIndex,1,1), sec_h, sec_b, longReinfNo, longReinfPhi, M_Rd, roMinCondition, shearPhi, shearSpac, shearLegs, V_Rd_it, sCondition, Fz_Ed];
    beamsMid(end+1,:) = [DataDesignMax(barIndex,1,1), sec_h, sec_b, longReinfNo, longReinfPhi, M_Rd, 0, shearPhiMid, shearSpacMid, shearLegsMid, V_RdMid, Fz_Ed];
    
    waitbar(size(beams,1) / length(beamDesiOrd),loading,'Beams progress','Name', 'DC1: Step 2 of 4');
end
save([folder '\DC1beamsIt1.mat'],   'beams')
save([folder '\DC1beamsIt1mid.mat'],'beamsMid')

clear barIndex beamsMid Fz_Ed i longReinfArea longReinfNo longReinfPhi M_Rd maxM_Ed roMinCondition sCondition ...
    sec_h sec_b shearLegs shearLegsMid shearPhi shearPhiMid shearSpac shearSpacMid V_Rd_it V_RdMid
%% STEP 3 - DESIGN COLUMNS
close(loading); loading = waitbar(0,'Initializing columns','Name', 'DC1: Step 3 of 4'); pause(1);

noStories = max(stories(:,1)); count = 0; columns = []; columnsMid = [];
tic
for i = 1 : size(barsOfColumns,1)
    
    barNames = []; % to append in the end
    auxColumns = [];
    
    for j = noStories : -1 : 1
%         display(j)
        barName = barsOfColumns(i,j); barNames = [barNames; barName];
        barIndex = find(DataDesignMax(:,1,1) == barName);
        
        try minWidth = minWidFind(barName, element, beams);
            catch minWidth = .2;
        end
        try minWidth = max(max(minWidth, auxColumns(:,1)));
            catch minWidth = .2;
        end
        
        % longitudinal reinforcement
        mAux1 = [];
        for k = allCasesIdx
            N_axial = DataDesignMax(barIndex, 2, k);
            My_h    = DataDesignMax(barIndex, 5, k);
            Mz_b    = DataDesignMax(barIndex, 6, k);
            V_Ed    = 'skip';
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(N_axial, My_h, Mz_b, V_Ed, minWidth);
            mAux1 = [mAux1; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, 0]];
            
            N_axial = DataDesignMin(barIndex, 2, k);
            My_h    = DataDesignMin(barIndex, 5, k);
            Mz_b    = DataDesignMin(barIndex, 6, k);
            V_Ed    = 'skip';
            [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(N_axial, My_h, Mz_b, V_Ed, minWidth);
            mAux1 = [mAux1; [sec_h, sec_b, noRebar, phiRebar, areaRebar, reinfPercFin, M_Rd, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, 0]];
        end
        % flex dims and long
        mAux1 = sortrows(mAux1,[5 1],{'descend' 'descend'});
        
        % shear columns reinforcement
        givenWidth = mAux1(1, 1); givenLong  = mAux1(1, [4, 3, 5]);
        
        V_Ed = max(max(abs(DataDesignMax(barIndex, 4, allCasesIdx)), abs(DataDesignMin(barIndex, 4, allCasesIdx))));
        
        N_axial = 'skip'; My_h = 'skip'; Mz_b = 'skip';
        [~, ~, ~, ~, ~, ~, ~, shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition] = DC1columnDesign(N_axial, My_h, Mz_b, V_Ed, givenWidth, givenLong);

        % finish design
        auxColumns(j, :)      = mAux1(1, :);
        auxColumns(j, [8:14]) = [shearReinfPhi, shearReinfSpac, shearReinfLoops, shearReinfArea, V_Rd, sCondition, V_Ed];
        toc
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    
    columns = [columns; [barNames, auxColumns]];
    %     columnsMid = [columnsMid; [barNames, auxColumns(:,[1:7]), auxColumnsMid]];
    
%     display(columns)%%%%%%%%%%%%%%%%%%%%%%%%%
    
    count = count + 1;
    loading = waitbar(count / size(barsOfColumns,1),loading,'Columns progress','Name', 'DC1: Step 3 of 4');
end
save([folder '\DC1columnsIt1.mat'],   'columns')
save([folder '\DC1columnsIt1mid.mat'],'columnsMid')
%% STEP 4 - EXPORT RESULTS TO CSV FOR SEISMOSTRUCT
close(loading); loading = waitbar(0,'Updating and writing to Seismo','Name', 'DC1: Step 4 of 4'); pause(1);
toSeismo(columns(:,[1:5, 9:11]), beams(:,[1:5, 8:10]), nodes, element, stories, DataDesignMin, DataDesignMax, seismicVerticalLoadCaseIdx, fck, fyk, cover, folder);
loading = waitbar(1, loading,'Updating and writing to Seismo','Name', 'DC1: Step 4 of 4');
%% FINISH
time = toc; save([folder '\time.mat'],'time');
minutes = floor(time/60);
seconds = floor(time - minutes*60);
disp(['Finished in ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
close(loading);