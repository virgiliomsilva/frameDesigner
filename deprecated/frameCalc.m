
function frameCalc(nFloor,nBay,dsgAccel,fHeight,spanL)
    %% design of structures
    % this is a self-contained script only requiring OpenSees to work

    % OpenSees' global axis for 3D models 
    %       | Y
    %       |
    %       |
    %       |_______ X
    %      /
    %     /
    %    /Z
    %% initial check 
    if exist('OpenSees.exe','file')==0;
    error('Error 404: OpenSees.exe not found.');
    end
    %%
    if nargin==3;
         % in the absence of user defined floor height or span
         % the script will randonly generate these values from the distributions
         % in (Silva et al 2013)
        randAux=normrnd([2.88,log(4.37)-1/2*(log(1+0.11^2))],[0.07*2.88,sqrt((log(1+0.11^2)))]);
        fHeight1=max([2.7 round(randAux(1)*10)/10]); fHeight=min([fHeight1 3.2]);
        spanL=round(exp(randAux(2))*10)/10;
    end
    nBayX=nBay;nBayZ=nBay;
    qk=2.0+0.8; % kN/m^2 (live loads table 6.2 + 6.3.1.2(8)P (EN1991-1))
    qkTop=0.4; % kN/m^2 (live load for top floor table 6.10 (EN1991-1))
    gk=25*0.25; % 0.25m RC slab (kN/m^2)
    fck=25; % MPa
    cover=0.03; % [m]
    fy=500; % Mpa
    fctm=0.3*fck^(2/3);fcm=fck+8;
    Ec=5000*sqrt(fcm);Gc=Ec/(2*(1+0.2));
    epsc0=2/1000;epscU=5/1000;
    clear randAux    
    %% --------------------------(end 1st stage)--------------------------
    %% write auxiliar tcl files
    %----
    if exist('staticAnalysis.tcl','file')==0;
        fid=fopen('staticAnalysis.tcl','wt+');
        fprintf(fid,'set nSteps 21');fprintf(fid,'\n');
        fprintf(fid,'for {set iStep 1} {$iStep<$nSteps} {incr iStep 1} {');
        fprintf(fid,'\n');fprintf(fid,'system BandGeneral');
        fprintf(fid,'\n'); fprintf(fid,'constraints Transformation');
        fprintf(fid,'\n');fprintf(fid,'numberer RCM');fprintf(fid,'\n');
        fprintf(fid,'test NormDispIncr 1.0e-5  10 3');fprintf(fid,'\n');
        fprintf(fid,'algorithm Newton');fprintf(fid,'\n');fprintf(fid,'integrator LoadControl 0.05');
        fprintf(fid,'\n');fprintf(fid,'analysis Static');fprintf(fid,'\n');
        fprintf(fid,'set ok [analyze 1]');fprintf(fid,'\n');
        fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
        fprintf(fid,'\n'); fprintf(fid,'algorithm Newton -initial');fprintf(fid,'\n');
        fprintf(fid,'set ok [analyze 1]}');fprintf(fid,'\n');
        fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
        fprintf(fid,'\n'); fprintf(fid,'algorithm ModifiedNewton -initial');fprintf(fid,'\n');
        fprintf(fid,'set ok [analyze 1]}}');fprintf(fid,'\n');
        fprintf(fid,'loadConst -time 0.0');
        fclose all;
    end
    %----
    %% psi values for load combinations 
    psi0=[0.6 0.7]; %[psi0Wind psi0LiveLoad]

    psi2=[0 0.3]; %[psi2Wind psi2LiveLoad]

    %% compute nodes' coordinates, tag fixed nodes and compute masses
    % for final 3D frame
    nCoord3D=zeros((nBayZ+1)*((nBayX+1)*(nFloor+1)),4); %[nTag xCoord yCoord zCoord]
    nNodes3D=0;ndFlag=(zeros((nBayZ+1)*((nBayX+1)*(nFloor+1)),1));
    fixN3D=(zeros((nBayZ+1)*((nBayX+1)*(nFloor+1)),1));
    for iBayZ=1:(nBayZ+1);
        for iBayX=1:(nBayX+1);
            for iFloor=1:(nFloor+1);
                nNodes3D=nNodes3D+1;
                if iFloor==1;
                    fixN3D(nNodes3D,1)=1;
                end
                xCoord=(iBayX-1)*spanL;
                yCoord=(iFloor-1)*fHeight;
                zCoord=(iBayZ-1)*spanL;
                ndFlag(nNodes3D,1)=str2num([num2str(iBayZ) num2str(iBayX) num2str(iFloor)]);
                nCoord3D(nNodes3D,1)=nNodes3D;
                nCoord3D(nNodes3D,2)=xCoord;
                nCoord3D(nNodes3D,3)=yCoord;
                nCoord3D(nNodes3D,4)=zCoord;
            end
        end
    end

    % for auxiliar 2D frames (used for initial design)
    nCoord2D=zeros((nBay+1)*(nFloor+1),3);
    fixN2D=zeros((nBay+1)*(nFloor+1),1);
    nNodes2D=0;
    for iBay=1:nBay+1;
        for iFloor=1:nFloor+1;
            nNodes2D=nNodes2D+1;
            xCoord=(iBay-1)*spanL;yCoord=(iFloor-1)*fHeight;
            nCoord2D(nNodes2D,1)=nNodes2D;
            nCoord2D(nNodes2D,2)=xCoord;
            nCoord2D(nNodes2D,3)=yCoord;
            if iFloor==1;
                fixN2D(nNodes2D,1)=1;
            end
        end
    end

    % compute masses  
    iflArea(1,1)=(spanL/2)^2;iflArea(2,1)=(spanL/2)*spanL;iflArea(3,1)=spanL^2;
    %[corner middleFrame centreStructure]
    mass3D=zeros(nNodes3D,1); count=0;
    for iBayZ=1:(nBayZ+1);
        for iBayX=1:(nBayX+1);
            for iFloor=1:(nFloor+1);
                count=count+1;
                ped=gk+(psi2(2)*qk);
                if iFloor==nFloor+1;
                    ped=gk+(psi2(2)*qkTop);
                end
                if iFloor>1;
                    if iBayZ==1||iBayZ==(nBayZ+1);
                        if iBayX==1||iBayX==(nBayX+1);
                            mass3D(count,1)=1/9.81*(ped*iflArea(1,1));
                        else
                            mass3D(count,1)=1/9.81*(ped*iflArea(2,1));
                        end
                    else
                        if iBayX==1||iBayX==(nBayX+1);
                            mass3D(count,1)=1/9.81*(ped*iflArea(2,1));
                        else
                            mass3D(count,1)=1/9.81*(ped*iflArea(3,1));
                        end
                    end
                end                    
            end
        end    
    end
    massTot3D=sum(mass3D(:,1));

    mass2D=zeros(nNodes2D,2);% mass for [extFrame intFrame]
    count=0;
    for iBay=1:(nBay+1);
        for iFloor=1:nFloor+1;
            ped=gk+psi2(2)*qk;
            count=count+1;
            if iFloor==nFloor+1;
                ped=gk+psi2(2)*qkTop;
            end
            if iFloor>1;
                if iBay==1||iBay==(nBay+1);
                    mass2D(count,1)=1/9.81*(ped*iflArea(1,1));
                    mass2D(count,2)=1/9.81*(ped*iflArea(2,1));

                else
                    mass2D(count,1)=1/9.81*(ped*iflArea(2,1));
                    mass2D(count,2)=1/9.81*(ped*iflArea(3,1));
                end
            end           
        end
    end
    massTot2D=[sum(mass2D(:,1)) sum(mass2D(:,2))];

    clear xCoord yCoord zCoord iBay iBayX iBayZ iFloor count

    %% Preliminary design 

    %------------ (vertical loads + deformation considerations) ---------------

    % compute column sections (1st Iteration-vertical loads only)

    secCol=zeros(size(iflArea,1),6);% matrix to store column properties [B H n phi phiT sT]
    for iSec=1:size(iflArea,1);    
        fcd=(fck*1000)/1.5;
        nsd=(1.35*gk+1.5*qk)*iflArea(iSec,1)*1.1*(nFloor-1)+...
            (1.35*gk+1.5*qkTop)*iflArea(iSec,1)*1.1;% axial load for predesign (kN)
        %nSd=fcd*Ac+fy2/1000*As
        ac=nsd/(fcd+400e3*0.01);% As approx 0.01Ac
        % compute section dimentions
        % round to the nearest multiple of 0.05m and min of 0.25m
        secDim=max([0.25 round(sqrt(ac)/0.05)*0.05]);
        secCol(iSec,1)=secDim;secCol(iSec,2)=secDim;
    end
    clear nsd ac as asR phi nB nRebar secDim iSec iPhi auxLogical index indexV

    % beams sections (1st Iteration-vertical loads only)
    % ACI method for pre-design
    % cross reinforced slab->A=b*h/2
    % ---
    % compute positiveM
    mBeamPos=zeros(nBayZ+1,nBayX);
    for iBayX=1:nBayX;
        for iBayZ=1:nBayZ+1;
            if iBayZ==1||iBayZ==nBayZ+1; %number of frames in Z
                iflAreaB=(spanL*(spanL/2))/2;
            else
                iflAreaB=((spanL*(spanL/2))/2)*2;
            end
                psd=((1.5*gk+1.5*qk)*iflAreaB)/spanL;        
            if iBayX==1||iBayX==nBayX;
                mBeamPos(iBayZ,iBayX)=(psd*spanL^2)/14;
            else
                mBeamPos(iBayZ,iBayX)=(psd*spanL^2)/16;
            end
        end
    end
    % compute negativeM
    mBeamNeg=zeros(nBayZ+1,nBayX);
    for iBayX=1:nBayX;
        for iBayZ=1:nBayZ+1;
            if iBayZ==1||iBayZ==nBayZ+1; %number of frames in Z
                iflAreaB=(spanL*(spanL/2))/2;
            else
                iflAreaB=((spanL*(spanL/2))/2)*2;
            end
                psd=((1.5*gk+1.5*qk)*iflAreaB)/spanL;
            if iBayX==1||iBayX==nBayX;
                if nBay>2;
                    mBeamNeg(iBayZ,iBayX)=-(psd*spanL^2)/10;
                else
                    mBeamNeg(iBayZ,iBayX)=-(psd*spanL^2)/9;
                end
            else
                mBeamNeg(iBayZ,iBayX)=-(psd*spanL^2)/11;
            end        
        end
    end
    clear iBayX iBayZ psd

    % compute 1st dimensions for beams (b approx 0.4*h)
    secBeam=zeros(2,12); %[b H nRebarUPM- PhiUPM- nRebarDownM- PhiDownM- nRebarUPM+ PhiUPM+ nRebarDownM+ PhiDownM+ phiT sT] 
    for iBay=1:2;
        msd=max(abs([mBeamPos(iBay,:) mBeamNeg(iBay,:)]));
        d=(msd/(0.1*(fck*1000)/1.5))^(1/3);
        h=max([ceil((d+cover)/0.05)*0.05 ceil((spanL/12)/0.05)*0.05]);
        b=max([0.2 ceil((0.4*h)/0.05)*0.05]);

        secBeam(iBay,1)=b;
        secBeam(iBay,2)=h;      
    end
    clear d b h iBay
    %% 
    controlVal=0;nIter=0;
    while controlVal==0 && nIter<=10;
        %%  compute period of vibration and eigenvector with elastic 2D models 
        if nBay>1;
            frame={'Ext','Int'};% exterior and interior frames
        else
            frame={'Ext'};
        end
        T1Vec=zeros(1,size(frame,2));phi1Vec=cell(1,size(frame,2));
        for iFrame=1:size(frame,2);
            fileName=char(strcat('frame2D',frame(iFrame),'.tcl'));
            folderName=char(strcat('OutFrame2D',frame(iFrame)));
            if exist(folderName,'dir')==0;
                outFolder=mkdir(folderName);
            end
            fid=fopen(fileName,'wt+');
            fprintf(fid,'wipe all');fprintf(fid,'\n');
            fprintf(fid,'set dataDir %s',folderName);fprintf(fid,'\n');
            fprintf(fid,'logFile $dataDir/logfile.txt');fprintf(fid,'\n');
            fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3');
            fprintf(fid,'\n');
            % create nodes
            for iNode=1:nNodes2D;
                fprintf(fid,'node %g %g %g \t',nCoord2D(iNode,1),nCoord2D(iNode,2),...
                    nCoord2D(iNode,3));fprintf(fid,'\n');
                if fixN2D(iNode,1)==1;
                    fprintf(fid,'fix %g %g %g %g \t',nCoord2D(iNode,1),1,1,1);
                    fprintf(fid,'\n'); 
                else
                fprintf(fid,'mass %g %g %g %g \t',nCoord2D(iNode,1),...
                    mass2D(iNode,iFrame),mass2D(iNode,iFrame),1e-9);
                fprintf(fid,'\n');
                end
            end
            % write geomtransf
            fprintf(fid,'geomTransf Corotational 1 \t');
            fprintf(fid,'\n');
            fprintf(fid,'geomTransf Linear 2 \t');
            fprintf(fid,'\n');
            % create elastic elements
            % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
            eleTag=0;elem2D=[]; % matrix to store element connectivity [eleTag nodeI nodeJ]
            % columns
            for iBay=1:nBay+1;
                for iFloor=1:nFloor;
                    eleTag=eleTag+1;
                    xCoordI=(iBay-1)*spanL;yCoordI=(iFloor-1)*fHeight;
                    xCoordJ=(iBay-1)*spanL;yCoordJ=(iFloor)*fHeight;
                    auxMat1=nCoord2D(nCoord2D(:,2)==xCoordI,:);
                    auxMat2=nCoord2D(nCoord2D(:,2)==xCoordJ,:);
                    nodeI=auxMat1((find(auxMat1(:,3)==yCoordI)),1);
                    nodeJ=auxMat2((find(auxMat2(:,3)==yCoordJ)),1);
                    elem2D=[elem2D; [eleTag nodeI nodeJ]];
                    if iFrame==1;
                        if iBay==1||iBay==(nBay+1);
                            bCol=secCol(1,1);hCol=secCol(1,2);
                        else
                            bCol=secCol(2,1);hCol=secCol(2,2);
                        end
                    else
                        if iBay==1||iBay==(nBay+1);
                            bCol=secCol(2,1);hCol=secCol(2,2);
                        else
                            bCol=secCol(3,1);hCol=secCol(3,2);
                        end            
                    end
                    IzCol=(bCol*hCol^3)/12;ACol=bCol*hCol;            
                    fprintf(fid,'element elasticBeamColumn %g %g %g %g %g %g %g \t',...
                        eleTag,nodeI,nodeJ,ACol,Ec*1000,IzCol,1);fprintf(fid,'\n');          
                end
            end
            endColTag=eleTag;% tag to mark the end of the column elements
            % beams
            for iFloor=1:nFloor;
                for iBay=1:nBay;
                    eleTag=eleTag+1;
                    xCoordI=(iBay-1)*spanL;yCoordI=(iFloor)*fHeight;
                    xCoordJ=(iBay)*spanL;yCoordJ=(iFloor)*fHeight;
                    auxMat1=nCoord2D(nCoord2D(:,2)==xCoordI,:);
                    auxMat2=nCoord2D(nCoord2D(:,2)==xCoordJ,:);
                    nodeI=auxMat1((find(auxMat1(:,3)==yCoordI)),1);
                    nodeJ=auxMat2((find(auxMat2(:,3)==yCoordJ)),1);
                    elem2D=[elem2D; [eleTag nodeI nodeJ]];
                     if iFrame==1;
                            bBeam=secBeam(1,1);hBeam=secBeam(1,2);
                     else
                            bBeam=secBeam(2,1);hBeam=secBeam(2,2);
                     end      
                    IzBeam=(bBeam*hBeam^3)/12;ABeam=bBeam*hBeam;            
                    fprintf(fid,'element elasticBeamColumn %g %g %g %g %g %g %g \t',...
                        eleTag,nodeI,nodeJ,ABeam,Ec*1000,IzBeam,2);fprintf(fid,'\n');
                end
            end    
            %compute fundamental period and eigenvector
            fprintf(fid,'set pi %g \t',pi);fprintf(fid,'\n');
            fprintf(fid,'set lambda [eigen  1]');fprintf(fid,'\n');
            fprintf(fid,'set T [expr 2.0*$pi/pow([lindex $lambda 0],0.5)]');
            fprintf(fid,'\n');fprintf(fid,'set outfile [open "$dataDir/T1.txt" w]');
            fprintf(fid,'\n');fprintf(fid,'puts $outfile $T');fprintf(fid,'\n');
            fprintf(fid,'close $outfile');fprintf(fid,'\n');
            ctrlNode1=min(nCoord2D(nCoord2D(:,3)==fHeight,1));
            ctrlNode2=min(nCoord2D(nCoord2D(:,3)==(fHeight*nFloor),1));
            fprintf(fid,'recorder Node -file $dataDir/phi1.txt -nodeRange %g %g -dof 1 "eigen1" \t',...
                ctrlNode1,ctrlNode2);fprintf(fid,'\n');

            fprintf(fid,'source staticAnalysis.tcl');fprintf(fid,'\n');
            % close file 
            fclose all;    
            % 1st run of opensees
            eval(['!' 'OpenSees' ' ' fileName]);

            % import and store period and eigenvector in Matlab    
            T1Vec(1,iFrame)=importdata(char(strcat(folderName,'/T1.txt')),'\t');
            auxMat1=importdata(char(strcat(folderName,'/phi1.txt')));
            phi1Vec{1,iFrame}=auxMat1(end,:)';      
        end
        clear iFrame fileName folderName fid eleTag iNode xCoordI xCoordJ yCoordI...
            yCoordJ auxMat auxMat1 iBay iFloor bBeam hBeam bCol hCol...
            bLoad iflAreaB  Abeam Acol IzBeam IzCol

        nElem=size(elem2D,1);

        %% compute wind loads with EC1 (EN1991-4)
        % --- partI (constant values)---
        % (assumed terrain CatII)
        vb=25; %m/s
        z0=0.05;zMin=2;
        zVec=(0:fHeight:(nFloor*fHeight))';kr=0.19*(z0/0.05)^0.07;
        cr=zeros(size(zVec,1),1);
        for iLine=1:size(cr,1);
            if zVec(iLine,1)<=zMin;
                cr(iLine,1)=kr*log(zMin/z0);
            else
                cr(iLine,1)=kr*log(zVec(iLine,1)/z0);
            end    
        end
        vm=cr*vb;sigmaV=kr*vb;
        Iv=sigmaV./vm;qp=(1+7*Iv).*((1/2)*1.25*vm.^2);
        zs=0.6*(nFloor*fHeight);bSquare=1.0;
        % turbulence scale
        alpha=0.67+0.05*log(z0);Lt=300;zt=200;
        Lz=zeros(size(zVec,1),1);
        for iLine=1:size(zVec,1);
            if zVec(iLine,1)>=zMin;
                Lz(iLine)=Lt*(zVec(iLine,1)/zt)^alpha;
            else
                Lz(iLine)=Lt*(zMin/zt)^alpha;
            end    
        end
        clear iLine
        cf0Mat=[0.1 0.2 0.6 0.7 1 2 5 10 20 50;
            2 2 2.35 2.4 2.1 1.65 1 0.9 0.9 0.9]';

        % ---compute wind forces (Part II)---
        h=nFloor*fHeight; b=nBay*spanL; d=b;
        % compute pessure coefficients for vertical walls (table 7.1 EN1991-4)
        if h/d<=0.25;
            cpe10D=0.7;cpe10E=-0.3;
        elseif h/d<=1;
            cpe10D=0.7+(0.1/0.75)*(h/d-0.25);
            cpe10E=(-1*0.3)-(0.2/0.75)*(h/d-0.25);
        elseif h/d<=5
            cpe10D=0.8;
            cpe10E=(-1*0.5)-(0.2/4)*(h/d-0.25);
        else
            cpe10D=0.8;cpe10E=-0.7;
        end
        latWindF=cell(1,size(frame,2));%cell to store lateral wind forces
        topWindF=cell(1,size(frame,2));%cell to store roof wind forces
        for iFrame=1:size(frame,2);
            T1=T1Vec(1,iFrame);
            f1=1/T1;
            fl=f1*(Lz./vm);Sl=(6.8*fl)./((1+10.2*fl).^(5/3));
            Lzs=Lt*(zs/zt)^alpha;vmzs=(kr*log(zs/z0))*vb;
            flzs=f1*(Lzs/vmzs);Slzs=(6.8*flzs)/((1+10.2*flzs).^(5/3));
            etaH=4.6*(nFloor*fHeight)/Lzs*flzs;
            etaB=4.6*(nBay*spanL)/Lzs*flzs;
            Rh=1/etaH-1/(2*etaH^2)*(1-exp(-2*etaH));
            Rb=1/etaB-1/(2*etaB^2)*(1-exp(-2*etaB));
            aux=0; aux2=0; phi1=phi1Vec{1,iFrame}(:,1);
            for iFloor=1:nFloor;
                index=find(nCoord2D(:,3)==(iFloor*fHeight));

                aux=aux+(sum(mass2D(index,iFrame)*1000)*(phi1(iFloor,1))^2);
                aux2=aux2+(phi1(iFloor,1))^2;        
            end
            me=aux/aux2;
            deltaA=2.1*1.25*(spanL*nBay)*vmzs/(2*f1*me);deltaS=0.10;
            rSquare=pi^2/(2*deltaA+deltaS)*Slzs*Rh*Rb;
            nu=max([0.08 f1*sqrt(rSquare/(bSquare+rSquare))]);
            kp=max([3 sqrt(2*log(nu*600))+0.6/(sqrt(2*log(nu*600)))]);
            Ivzs=sigmaV./vmzs;
            % compute cscd
            cscd=(1+2*kp*Ivzs*sqrt(bSquare+rSquare))/(1+7*Ivzs);
            % compute wind forces
            if iFrame==1;
                aux1=spanL/2;
            else
                aux1=spanL;
            end
            % --- compute lateral loads ---
            windNodalLoad=zeros(nFloor+1,2); % [forceZoneD forceZoneE] in Newton
            for iFloor=1:nFloor+1;
                cHeight=(iFloor-1)*fHeight;%current height of structural node
                heightUp=cHeight+fHeight/2;
                heightDown=max([0 cHeight-fHeight/2]); 
                if h<=b;
                    ze=h;
                    vmze=kr*log(ze/z0)*vb;
                    Ivze=sigmaV./vmze;
                    qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                    windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                    windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                elseif b<h && h<=(2*b);
                    if heightUp<=b && heightDown<b;
                        ze=b;
                        vmze=kr*log(ze/z0)*vb;
                        Ivze=sigmaV./vmze;
                        qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                    elseif heightUp>b && heightDown<=b;
                        ze1=h;
                        vmze1=kr*log(ze1/z0)*vb;
                        Ivze1=sigmaV./vmze1;
                        qpze1=(1+7*Ivze)*(1/2*1.25*vmze1^2);
                        ze2=b;
                        vmze2=kr*log(ze2/z0)*vb;
                        Ivze2=sigmaV./vmze2;
                        qpze2=(1+7*Ivze2)*(1/2*1.25*vmze2^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze1*cpe10D*(aux1*(heightUp-b))+...
                            cscd*2.1*qpze2*cpe10D*(aux1*(b-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze1*cpe10E*(aux1*(heightUp-b))+...
                            cscd*2.1*qpze2*cpe10E*(aux1*(b-heightDown));
                    elseif heightUp>b && heightDown>b;
                        ze=h;
                        vmze=kr*log(ze/z0)*vb;
                        Ivze=sigmaV./vmze;
                        qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                    end            
                elseif h>(2*b);
                     if heightUp<=b && heightDown<=b;
                        ze=b;
                        vmze=kr*log(ze/z0)*vb;
                        Ivze=sigmaV./vmze;
                        qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                     elseif heightUp<=(h-b) && heightDown<=b;
                        ze1=abs(heightUp+(h-b))/2;
                        vmze1=kr*log(ze1/z0)*vb;
                        Ivze1=sigmaV./vmze1;
                        qpze1=(1+7*Ivze1)*(1/2*1.25*vmze1^2);
                        ze2=b;
                        vmze2=kr*log(ze2/z0)*vb;
                        Ivze2=sigmaV./vmze2;
                        qpze2=(1+7*Ivze2)*(1/2*1.25*vmze2^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze1*cpe10D*(aux1*(heightUp-b))+...
                            cscd*2.1*qpze2*cpe10D*(aux1*(b-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze1*cpe10E*(aux1*(heightUp-b))+...
                            cscd*2.1*qpze2*cpe10E*(aux1*(b-heightDown));
                     elseif heightUp<=(h-b) && heightDown>b;
                        ze=cHeight;
                        vmze=kr*log(ze/z0)*vb;
                        Ivze=sigmaV./vmze;
                        qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                     elseif heightUp>(h-b) && heightDown<=(h-b);
                        ze1=h;
                        vmze1=kr*log(ze1/z0)*vb;
                        Ivze1=sigmaV./vmze1;
                        qpze1=(1+7*Ivze1)*(1/2*1.25*vmze1^2);
                        ze2=abs(heightDown+(h-b))/2;
                        vmze2=kr*log(ze2/z0)*vb;
                        Ivze2=sigmaV./vmze2;
                        qpze2=(1+7*Ivze2)*(1/2*1.25*vmze2^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze1*cpe10D*(aux1*(heightUp-(h-b)))+...
                            cscd*2.1*qpze2*cpe10D*(aux1*((h-b)-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze1*cpe10E*(aux1*(heightUp-(h-b)))+...
                            cscd*2.1*qpze2*cpe10E*(aux1*((h-b)-heightDown));
                     elseif heightUp>(h-b) && heightDown>(h-b);
                        ze=h;
                        vmze=kr*log(ze/z0)*vb;
                        Ivze=sigmaV./vmze;
                        qpze=(1+7*Ivze)*(1/2*1.25*vmze^2);
                        windNodalLoad(iFloor,1)=cscd*2.1*qpze*cpe10D*(aux1*(heightUp-heightDown));
                        windNodalLoad(iFloor,2)=cscd*2.1*qpze*cpe10E*(aux1*(heightUp-heightDown));
                     end       
                end
            end
            latWindF{1,iFrame}=windNodalLoad/1000; % store lateral wind forces in kN
            % --- compute loads for roof ---
            ze=h; vmze=kr*log(ze/z0)*vb; Ivze=sigmaV./vmze;
            qpze=(1+7*Ivze)*(1/2*1.25*vmze^2); e=min([b 2*h]);
            windTopLoad=zeros(2,nBay);
            for iBay=1:nBay;
                left=(iBay-1)*spanL; right=iBay*spanL;
                if left<=e/10 && right<=e/10;
                    if iFrame==1;                
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-1.8*(aux1*(right-left));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-1.8*(aux1*(right-left));
                    else
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-1.2*(aux1*(right-left));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-1.2*(aux1*(right-left));
                    end
                elseif left<=e/10 && right<=e/2;
                    if iFrame==1;                
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-1.8*(aux1*(e/10-left))+...
                            cscd*2.1*qpze*-0.7*(aux1*(right-e/10));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-1.8*(aux1*(e/10-left))+...
                            cscd*2.1*qpze*-0.7*(aux1*(right-e/10));
                    else
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-1.2*(aux1*(e/10-left))+...
                            cscd*2.1*qpze*-0.7*(aux1*(right-e/10));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-1.2*(aux1*(e/10-left))+...
                            cscd*2.1*qpze*-0.7*(aux1*(right-e/10));
                    end
                elseif left<e/2 && right<=e/2;
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-0.7*(aux1*(right-left));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-0.7*(aux1*(right-left));
                elseif left<e/2 && right>e/2;
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-0.7*(aux1*(e/2-left))+...
                            cscd*2.1*qpze*-0.2*(aux1*(right-e/2));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*-0.7*(aux1*(e/2-left))+...
                            cscd*2.1*qpze*+0.2*(aux1*(right-e/2));
                else
                        windTopLoad(1,iBay)=cscd*2.1*qpze*-0.2*(aux1*(right-left));
                        windTopLoad(2,iBay)=cscd*2.1*qpze*+0.2*(aux1*(right-left));           
                end
                topWindF{1,iFrame}=((windTopLoad/1000)*-1)/spanL;% store wind forces in kN/m and in OS global axis 

            end
        end
        clear iFloor aux aux1 aux2 h b d e windNodalLoad windTopLoad left right...
            heightUp heightDown iFrame iBay index
        %% Compute seismic loads
        if dsgAccel>0;% only enters if design accel greater than 0
            latEQF=cell(1,size(frame,2));
            S=1; TB=[0.15 0.05]; TC=[0.4 0.25]; TD=[2 1.2];
            % comb=Gk+sum(psi2Qk)EN1990-1
            for iFrame=1:size(frame,2);        
            % compute spectral acceleration (rock site assumed)
                T1=T1Vec(1,iFrame);auxVec=zeros(1,2);
                for iS=1:2;
                    if T1>=0 && T1<=TB(iS);
                        auxVec(1,iS)=dsgAccel*S*(1+T1/TB(iS)*(1*2.5-1));
                    elseif T1<=TC(iS);
                         auxVec(1,iS)=dsgAccel*S*1*2.5;
                    elseif T1<=TD(iS);
                        auxVec(1,iS)=dsgAccel*S*1*2.5*(TC(iS)/T1);
                    elseif T1<=4;
                        auxVec(1,iS)=dsgAccel*S*1*2.5*(TC(iS)*TD(iS)/T1^2);
                    end
                end
                [Saccel,index]=max(auxVec);Sdispl=Saccel*(T1/(2*pi))^2;
                %
                if nFloor==1;
                    alphaU=1.1;
                elseif nBay==1;
                    alphaU=1.2;
                else
                    alphaU=1.3;
                end
                %q=4.5*alphaU;
                q=3*alphaU;
                % design spectrum
                if T1>=0 && T1<=TB(index);
                    Sd=dsgAccel*S*(2/3+T1/TB(index)*(2.5*q-2/3));
                elseif T1<=TC(index);
                    Sd=dsgAccel*S*(2.5/q);
                elseif T1<=TD(index);
                    Sd=max([dsgAccel*S*1*2.5/q*(TC(index)/T1) 0.2*dsgAccel]);
                else
                   Sd=max([dsgAccel*S*1*2.5/q*(TC(index)*TD(index)/T1^2) 0.2*dsgAccel]);
                end
                % compute total base shear
                if T1<=2*TC(index)&&nFloor>2;
                    lambda=0.85;
                else
                    lambda=1;
                end
                Fb=Sd*massTot2D(1,iFrame)*lambda;
                massF=zeros(nFloor,1);
                for iFloor=1:nFloor;
                    index=find(nCoord2D(:,3)==fHeight*iFloor);
                    massF(iFloor,1)=sum(mass2D(index,iFrame));
                end
                phi1=phi1Vec{1,iFrame}(:,1);
                latEQF{1,iFrame}=Fb*(phi1.*massF)./...
                    (phi1'*massF);     
            end
        end
        %% load combinations

        if dsgAccel==0;
            nComb=(size(psi0,2)*size(topWindF,2))+3;% 2 from the cpe10 for Zone I on the roof * 2+3 combinations for vertical loads only
        else
            nComb=(size(psi0,2)*size(topWindF,2))+3+1; % 2 from the cpe10 for Zone I on the roof * 2+3 combinations for vertical loads only+1 for seismic
        end
        localF=cell(nComb,size(frame,2));% cell to store local forces
        N=cell(nComb,size(frame,2));V=cell(nComb,size(frame,2));
        M=cell(nComb,size(frame,2));
        bLoadCell=cell(nComb,size(frame,2));
        count=0;
        for ipsi0=1:size(psi0,2);
            for iWind=1:size(topWindF,2);
                count=count+1;
                for iFrame=1:size(frame,2);
                    folderName=char(strcat('OutFrame2D',frame(iFrame)));
                    mFile=char(strcat('frame2D',frame(iFrame),'.tcl')); % master file
                    fileName=char(strcat('frame2D',frame(iFrame),'Comb',num2str(count),'.tcl'));  
                    fid=fopen(mFile);fidd=fopen(fileName,'wt+');
                    while ~feof(fid) ; % reads the original till last line
                        tline=fgets(fid); %
                        if feof(fid)~=1;
                           fwrite(fidd,tline);
                        end
                    end
                    fclose all ;
                    fid=fopen(fileName,'at+');
                     % add vertical loads
                     fprintf(fid,'pattern Plain 1 Linear {');fprintf(fid,'\n');
                     eleTag=endColTag;
                     bLoad=zeros(nElem,1);
                     for iFloor=1:nFloor;
                         for iBay=1:nBay;
                             eleTag=eleTag+1;
                             if iFrame==1;
                                 iflAreaB=(spanL*(spanL/2))/2;
                             else
                                 iflAreaB=((spanL*(spanL/2))/2)*2;
                             end
                             if ipsi0==1;
                                 bLoad(eleTag,1)=(-1*(1.35*gk+1.5*qk)*iflAreaB)/spanL;
                                 if iFloor==nFloor;
                                     bLoad(eleTag,1)=(-1*(1.35*gk+1.5*qkTop)*iflAreaB)/spanL+...
                                         1.5*psi0(1,ipsi0)*topWindF{1,iFrame}(iWind,iBay);
                                 end
                             else 
                                 bLoad(eleTag,1)=(-1*(1.35*gk+1.5*qk*psi0(1,ipsi0))*iflAreaB)/spanL;
                                 if iFloor==nFloor;
                                     bLoad(eleTag,1)=(-1*(1.35*gk+1.5*qkTop*0)*iflAreaB)/spanL+...
                                         1.5*topWindF{1,iFrame}(iWind,iBay);
                                 end
                             end
                             fprintf(fid,'eleLoad -ele %g -type -beamUniform %g \t',...
                                 eleTag,bLoad(eleTag,1)); fprintf(fid,'\n');                 
                         end
                     end
                     bLoadCell{count,iFrame}=bLoad;
                     % add wind lateral loads
                     auxVec1=nCoord2D(nCoord2D(:,2)==0,1);
                     auxVec2=nCoord2D(nCoord2D(:,2)==(nBay*spanL),1);
                     for iLine=1:size(auxVec1,1);
                         if ipsi0==1;
                            fprintf(fid,'load %g %g %g %g \t',...
                                auxVec1(iLine,1),latWindF{1,iFrame}(iLine,1)*1.5*psi0(1,ipsi0),...
                                0,0);fprintf(fid,'\n');
                            fprintf(fid,'load %g %g %g %g \t',...
                                auxVec2(iLine,1),latWindF{1,iFrame}(iLine,2)*1.5*psi0(1,ipsi0),...
                                0,0);fprintf(fid,'\n');
                         else
                            fprintf(fid,'load %g %g %g %g \t',...
                                auxVec1(iLine,1),latWindF{1,iFrame}(iLine,1)*1.5,...
                                0,0);fprintf(fid,'\n');
                            fprintf(fid,'load %g %g %g %g \t',...
                                auxVec2(iLine,1),latWindF{1,iFrame}(iLine,2)*1.5,...
                                0,0);fprintf(fid,'\n');
                         end
                     end          
                    fprintf(fid,'}'); fprintf(fid,'\n'); 
                    % record element local forces
                    fprintf(fid,'recorder Element -file $dataDir/%s.txt -time -eleRange %g %g localForce \t',...
                        strcat('localForces',num2str(count)),elem2D(1,1),elem2D(end,1));fprintf(fid,'\n');
                    % run static analysis
                    fprintf(fid,'source staticAnalysis.tcl');fprintf(fid,'\n'); 
                    % model file closure
                    fclose all;
                    % run OpenSees
                    eval(['!' 'OpenSees' ' ' fileName]);
                    %store element local forces
                    auxMat=importdata(char(strcat(folderName,'/localForces',num2str(count),'.txt')));
                    localF{count,iFrame}=auxMat(end,:); % store only last step
                    Nlft=-1*(auxMat(end,2:6:end-5));Nrgt=auxMat(end,5:6:end-2);
                    Vlft=auxMat(end,3:6:end-4);Vrgt=-1*(auxMat(end,6:6:end-1));
                    Mlft=-1*(auxMat(end,4:6:end-3));Mrgt=(auxMat(end,7:6:end));
                    Mmid=((Mlft+Mrgt)./2)+(-bLoad'.*spanL^2)./8;
                    N{count,iFrame}=[Nlft' Nrgt'];V{count,iFrame}=[Vlft' Vrgt'];
                    M{count,iFrame}=[Mlft' Mmid' Mrgt'];
                    delete(fileName);
                end        
            end    
        end
        % (vertical loads only)
        for iComb=1:3;
            count=count+1;
            for iFrame=1:size(frame,2);
                folderName=char(strcat('OutFrame2D',frame(iFrame)));
                mFile=char(strcat('frame2D',frame(iFrame),'.tcl')); % master file
                fileName=char(strcat('frame2D',frame(iFrame),'Comb',num2str(count),'.tcl'));  
                fid=fopen(mFile);fidd=fopen(fileName,'wt+');
                while ~feof(fid) ; % reads the original till last line
                    tline=fgets(fid); %
                    if feof(fid)~=1;
                       fwrite(fidd,tline);
                    end
                end
                fclose all ;
                fid=fopen(fileName,'at+');
                % add vertical loads
                fprintf(fid,'pattern Plain 1 Linear {');fprintf(fid,'\n');
                eleTag=endColTag;
                bLoad=zeros(nElem,1);
                for iFloor=1:nFloor;
                     for iBay=1:nBay;
                         eleTag=eleTag+1;
                         if iFrame==1;
                             iflAreaB=(spanL*(spanL/2))/2;
                         else
                             iflAreaB=((spanL*(spanL/2))/2)*2;
                         end
                         if iFloor==nFloor
                             live=qkTop;
                         else 
                             live=qk;
                         end                     
                         if iComb==1;
                             bLoad(eleTag,1)=(-1*(1.35*gk+1.5*live)*iflAreaB)/spanL;
                         elseif iComb==2;
                             if mod(iFloor,2)~=0&&mod(iBay,2)~=0;
                                 bLoad(eleTag,1)=(-1*(1.35*gk+1.5*live)*iflAreaB)/spanL;
                             elseif mod(iFloor,2)==0&&mod(iBay,2)==0;
                                 bLoad(eleTag,1)=(-1*(1.35*gk+1.5*live)*iflAreaB)/spanL;
                             else
                                 bLoad(eleTag,1)=(-1*(1.35*gk)*iflAreaB)/spanL;
                             end
                         else 
                              if mod(iFloor,2)~=0&&mod(iBay,2)~=0;
                                 bLoad(eleTag,1)=(-1*(1.35*gk)*iflAreaB)/spanL;
                             elseif mod(iFloor,2)==0&&mod(iBay,2)==0;
                                 bLoad(eleTag,1)=(-1*(1.35*gk)*iflAreaB)/spanL;
                              else
                                 bLoad(eleTag,1)=(-1*(1.35*gk+1.5*live)*iflAreaB)/spanL;
                              end
                         end
                             fprintf(fid,'eleLoad -ele %g -type -beamUniform %g \t',...
                                 eleTag,bLoad(eleTag,1)); fprintf(fid,'\n'); 
                     end
                end
                bLoadCell{count,iFrame}=bLoad;
                fprintf(fid,'}'); fprintf(fid,'\n'); 
                % record element local forces
                fprintf(fid,'recorder Element -file $dataDir/%s.txt -time -eleRange %g %g localForce \t',...
                    strcat('localForces',num2str(count)),elem2D(1,1),elem2D(end,1));fprintf(fid,'\n');
                % run static analysis
                fprintf(fid,'source staticAnalysis.tcl');fprintf(fid,'\n'); 
                % model file closure
                fclose all;
                % run OpenSees
                eval(['!' 'OpenSees' ' ' fileName]);
                %store element local forces
                auxMat=importdata(char(strcat(folderName,'/localForces',num2str(count),'.txt')));
                localF{count,iFrame}=auxMat(end,:); % store only last step
                Nlft=-1*(auxMat(end,2:6:end-5));Nrgt=auxMat(end,5:6:end-2);
                Vlft=auxMat(end,3:6:end-4);Vrgt=-1*(auxMat(end,6:6:end-1));
                Mlft=-1*(auxMat(end,4:6:end-3));Mrgt=(auxMat(end,7:6:end));
                Mmid=((Mlft+Mrgt)./2)+(-bLoad'.*spanL^2)./8;
                N{count,iFrame}=[Nlft' Nrgt'];V{count,iFrame}=[Vlft' Vrgt'];
                M{count,iFrame}=[Mlft' Mmid' Mrgt'];
                delete(fileName);
            end
        end

        % seismic loads
        if dsgAccel>0;
            count=count+1;
            for iFrame=1:size(frame,2);
                folderName=char(strcat('OutFrame2D',frame(iFrame)));
                mFile=char(strcat('frame2D',frame(iFrame),'.tcl')); % master file
                fileName=char(strcat('frame2D',frame(iFrame),'Comb',num2str(count),'.tcl'));  
                fid=fopen(mFile);fidd=fopen(fileName,'wt+');
                while ~feof(fid) ; % reads the original till last line
                    tline=fgets(fid); %
                    if feof(fid)~=1;
                       fwrite(fidd,tline);
                    end
                end
                fclose all ;
                fid=fopen(fileName,'at+');
                % add vertical loads
                fprintf(fid,'pattern Plain 1 Linear {');fprintf(fid,'\n');
                eleTag=endColTag;
                bLoad=zeros(nElem,1);
                for iFloor=1:nFloor;
                     for iBay=1:nBay;
                         eleTag=eleTag+1;
                         if iFrame==1;
                             iflAreaB=(spanL*(spanL/2))/2;
                         else
                             iflAreaB=((spanL*(spanL/2))/2)*2;
                         end
                         if iFloor==nFloor
                             live=qkTop;
                         else 
                             live=qk;
                         end
                         bLoad(eleTag,1)=(-1*(gk+psi2(2)*live)*iflAreaB)/spanL;
                         fprintf(fid,'eleLoad -ele %g -type -beamUniform %g \t',...
                                 eleTag,bLoad(eleTag,1)); fprintf(fid,'\n');
                     end
                end
                bLoadCell{count,iFrame}=bLoad;
                % add earthquake lateral loads
                auxVec1=nCoord2D(nCoord2D(:,2)==0,1);auxVec1(1)=[];
                for iLine=1:size(auxVec1,1);
                     fprintf(fid,'load %g %g %g %g \t',...
                       auxVec1(iLine,1),latEQF{1,iFrame}(iLine,1),...
                          0,0);fprintf(fid,'\n');
                end          
                fprintf(fid,'}'); fprintf(fid,'\n'); 
                % record element local forces
                fprintf(fid,'recorder Element -file $dataDir/%s.txt -time -eleRange %g %g localForce \t',...
                    strcat('localForces',num2str(count)),elem2D(1,1),elem2D(end,1));fprintf(fid,'\n');
                % run static analysis
                fprintf(fid,'source staticAnalysis.tcl');fprintf(fid,'\n'); 
                % model file closure
                fclose all;
                % run OpenSees
                eval(['!' 'OpenSees' ' ' fileName]);
                %store element local forces
                auxMat=importdata(char(strcat(folderName,'/localForces',num2str(count),'.txt')));
                localF{count,iFrame}=auxMat(end,:); % store only last step
                Nlft=-1*(auxMat(end,2:6:end-5));Nrgt=auxMat(end,5:6:end-2);
                Vlft=auxMat(end,3:6:end-4);Vrgt=-1*(auxMat(end,6:6:end-1));
                Mlft=-1*(auxMat(end,4:6:end-3));Mrgt=(auxMat(end,7:6:end));
                Mmid=((Mlft+Mrgt)./2)+(-bLoad'.*spanL^2)./8;
                N{count,iFrame}=[Nlft' Nrgt'];V{count,iFrame}=[Vlft' Vrgt'];
                M{count,iFrame}=[Mlft' Mmid' Mrgt'];
                delete(fileName); 
            end
        end

        clear ipsi0 iFrame auxVec1 count auxVec2 auxMat iWind Nlft Nrgt Vlft...
            Vrgt Mlft Mmid Mrgt fileName mFile bLoad live
        %%  load envelop 
        NEnvMax=cell(1,size(frame,2));VEnvMax=cell(1,size(frame,2));
        MEnvMax=cell(1,size(frame,2));
        NEnvMin=cell(1,size(frame,2));VEnvMin=cell(1,size(frame,2));
        MEnvMin=cell(1,size(frame,2));
        VEnvAbsMax=cell(1,size(frame,2));

        for iFrame=1:size(frame,2);
            auxMatN1=zeros(size(elem2D,1),2);auxMatN2=zeros(size(elem2D,1),2);
            auxMatV1=zeros(size(elem2D,1),2);auxMatV2=zeros(size(elem2D,1),2);
            auxMatAbsV=zeros(size(elem2D,1),2);
            auxMatM1=zeros(size(elem2D,1),3);auxMatM2=zeros(size(elem2D,1),3);
            for iElem=1:size(elem2D,1);
                auxN=zeros(nComb,2);auxV=zeros(nComb,2);auxM=zeros(nComb,3);
                for iComb=1:nComb;
                    auxN(iComb,1)=N{iComb,iFrame}(iElem,1);
                    auxN(iComb,2)=N{iComb,iFrame}(iElem,2);
                    auxV(iComb,1)=V{iComb,iFrame}(iElem,1);
                    auxV(iComb,2)=V{iComb,iFrame}(iElem,2);
                    auxM(iComb,1)=M{iComb,iFrame}(iElem,1);
                    auxM(iComb,2)=M{iComb,iFrame}(iElem,2);
                    auxM(iComb,3)=M{iComb,iFrame}(iElem,3);
                end
                auxMatN1(iElem,1)=max(auxN(:,1));auxMatN1(iElem,2)=max(auxN(:,2));
                auxMatV1(iElem,1)=max(auxV(:,1));auxMatV1(iElem,2)=max(auxV(:,2));
                auxMatM1(iElem,1)=max(auxM(:,1));auxMatM1(iElem,2)=max(auxM(:,2));
                auxMatM1(iElem,3)=max(auxM(:,3));

                auxMatN2(iElem,1)=min(auxN(:,1));auxMatN2(iElem,2)=min(auxN(:,2));
                auxMatV2(iElem,1)=min(auxV(:,1));auxMatV2(iElem,2)=min(auxV(:,2));
                auxMatM2(iElem,1)=min(auxM(:,1));auxMatM2(iElem,2)=min(auxM(:,2));
                auxMatM2(iElem,3)=min(auxM(:,3));

                auxMatAbsV(iElem,1)=max(abs(auxV(:,1)));auxMatAbsV(iElem,2)=max(abs(auxV(:,2)));
            end
            NEnvMax{1,iFrame}=auxMatN1;NEnvMin{1,iFrame}=auxMatN2;
            VEnvMax{1,iFrame}=auxMatV1;VEnvMin{1,iFrame}=auxMatV2;
            MEnvMax{1,iFrame}=auxMatM1;MEnvMin{1,iFrame}=auxMatM2;
            VEnvAbsMax{1,iFrame}=auxMatAbsV;
        end

        clear auxN auxV auxM auxMatN1 auxMatV1 auxMatM1 auxMatN2 auxMatV2 auxMatM2...
            auxMat auxMat2
        %% design logitudinal rebar
        %  beam rebar
        calcAsBeam=cell(1,size(frame,2));
        for iFrame=1:size(frame,2);
            count=0;
            if fck<=50;
                k1=0.44;k2=1.25*(0.6+0.0014/3.5e-3);
                xOvD=0.25; % maximum x/d ratio
                deltaMax=max([0.7 k1+k2*xOvD]);
            else
                ecu2=mean([3.1 2.9 2.7 2.6 2.6])/1000;
                k3=0.54;k4=1.25*(0.6+0.0014/ecu2);
                xOvD=0.15; % maximum x/d ratio 
                deltaMax=max([0.7 k3+k4*xOvD]);
            end 
            auxMat=zeros(nElem-endColTag,4);% [AsUp AsDown miu]
            for iElem=(endColTag+1):nElem;
                count=count+1; hSec=secBeam(iFrame,2); bSec=secBeam(iFrame,1);
                dSec=(hSec-cover);
                auxVec=zeros(nComb,4);%[AsUpM- AsDownM- AsUpM+ AsDownM+]
                for iComb=1:nComb;
                    bLoad=bLoadCell{iComb,iFrame}(iElem,1);
                    MLft=M{iComb,iFrame}(iElem,1); 
                    MRgt=M{iComb,iFrame}(iElem,3);
                    MMid=M{iComb,iFrame}(iElem,2);
                    if iComb<8;
                       if sign(MLft)==sign(MRgt) && sign(MLft)==-1;
                           Mend=max(abs([MLft MRgt]))*deltaMax;
                           index=find(abs([MLft MRgt])<Mend);
                           if isempty(index)==1;
                               Mmid=-Mend+abs(bLoad*spanL^2/8);
                           else
                               [~,Maux]=find(abs([MLft MRgt])<Mend);

                               Mmid=((-1*(Mend)+Maux)/2)+abs(bLoad*spanL^2/8);
                           end
                       elseif sign(MLft)~=sign(MRgt) && sign(MLft)==1;
                           Mend=max([abs(MRgt)*deltaMax MLft]);
                           Mmid=(MLft-Mend)/2+abs(bLoad*spanL^2/8);
                       elseif sign(MLft)~=sign(MRgt) && sign(MLft)==-1
                           Mend=max([abs(MLft)*deltaMax MRgt]);
                           Mmid=(-1*Mend+MRgt)/2+abs(bLoad*spanL^2/8);                              

                       end
                    else
                        Mend=max(abs([MLft MRgt])); Mmid=abs(MMid);
                    end                     
                    asMin=max([0.226*fctm/(fy/1.15)*bSec*dSec 0.0013*bSec*dSec]);
                    tol=1;xPos=(cover-0.001);Med=Mend;
                    M1=0;
                    while tol>0.01&&xPos<hSec&&abs(M1)<Med;
                        xPos=xPos+0.001;
                        asCalc1=(0.8*xPos*bSec*(fck*1000/1.5))/((fy*1000/1.15)...
                            -200e6*(3.5e-3*(xPos-cover)/xPos));
                        M1=(0.8*xPos*bSec*(fck*1000/1.5))*(dSec-0.4*xPos)+...
                            asCalc1*200e6*(3.5e-3*(xPos-cover)/xPos)*(dSec-cover);
                        tol=abs((M1-Med)/Med);
                    end
                    tol=1;xPos=(cover-0.001);Med=Mmid;
                    M1=0;
                    while tol>0.01&&xPos<hSec&&abs(M1)<Med;
                        xPos=xPos+0.001;
                        asCalc2=(0.8*xPos*bSec*(fck*1000/1.5))/((fy*1000/1.15)...
                            -200e6*(3.5e-3*(xPos-cover)/xPos));
                        M1=(0.8*xPos*bSec*(fck*1000/1.5))*(dSec-0.4*xPos)+...
                            asCalc2*200e6*(3.5e-3*(xPos-cover)/xPos)*(dSec-cover);
                        tol=abs((M1-Med)/Med);
                    end
                    auxVec(iComb,1)=max([asCalc1 asMin/2]);
                    auxVec(iComb,2)=max([asCalc1 asMin/2]);
                    auxVec(iComb,3)=max([asCalc2 asMin/2]);    
                    auxVec(iComb,4)=max([asCalc2 asMin/2]);
                end            
            auxMat(count,1)=max(auxVec(:,1));auxMat(count,2)=max(auxVec(:,2));
            auxMat(count,3)=max(auxVec(:,3));auxMat(count,4)=max(auxVec(:,4));
            end
            calcAsBeam{1,iFrame}=auxMat;
        end
        MrdBeam=zeros(size(secBeam,1),1);MrdBeam2=zeros(size(secBeam,1),1);
        phiRebar=[0.012;0.016;0.02;0.025;0.032];% possible phi for rebar
         for iFrame=1:size(frame,2);
              as(1,1)=max(calcAsBeam{1,iFrame}(:,1));% AsUpM-
              as(2,1)=max(calcAsBeam{1,iFrame}(:,2));% AsDwnM-
              as(3,1)=max(calcAsBeam{1,iFrame}(:,3));% AsUpM+
              as(4,1)=max(calcAsBeam{1,iFrame}(:,4));% AsDwnM+
              auxVec1=(3:2:9);auxVec2=(4:2:10);
              for iLine=1:4;
                  if as(iLine,1)>0;
                      nRebar=max([2*ones((size(phiRebar,1)),1) ceil(as(iLine,1)./((pi.*phiRebar.^2)./4))],[],2);
                      asR1=nRebar.*(pi*phiRebar.^2./4);
                      [~,index1]=min((asR1-(ones([size(asR1,1),1])*as(iLine,1))));
                      if nRebar(index1)>=2
                        space=(secBeam(iFrame,1)-(2*cover))/(nRebar(index1)-1);
                        if space<phiRebar(index1)
                          index1=index1+1;
                        end
                      end
                      secBeam(iFrame,(auxVec1(iLine)))=nRebar(index1);
                      secBeam(iFrame,(auxVec2(iLine)))=phiRebar(index1);
                  end
              end 
         end

         for iFrame=1:size(frame,2);
            asBeam=secBeam(iFrame,3)*(pi*secBeam(iFrame,4)^2/4);
            bSec=secBeam(iFrame,1); hSec=secBeam(iFrame,2);
            dSec=hSec-cover;Fs1=asBeam*fy*1000/1.15;
            tol=1; xPos=cover-0.001;
            while tol>0.05 && xPos<hSec;
                xPos=xPos+0.001;        
                Fc=0.8*xPos*bSec*fck*1000/1.5;
                Fs2=asBeam*200e6*(3.5e-3*(xPos-cover)/xPos);
                tol=abs((Fs1-(Fc+Fs2))/Fs1);        
            end
            MrdBeam(iFrame,1)=Fc*(dSec-0.4*xPos)+Fs2*(dSec-cover);

            asBeam=secBeam(iFrame,7)*(pi*secBeam(iFrame,8)^2/4);
            bSec=secBeam(iFrame,1); hSec=secBeam(iFrame,2);
            dSec=hSec-cover;Fs1=asBeam*fy*1000/1.15;
            tol=1; xPos=cover-0.001;
            while tol>0.05 && xPos<hSec;
                xPos=xPos+0.001;        
                Fc=0.8*xPos*bSec*fck*1000/1.5;
                Fs2=asBeam*200e6*(3.5e-3*(xPos-cover)/xPos);
                tol=abs((Fs1-(Fc+Fs2))/Fs1);        
            end
            MrdBeam2(iFrame,1)=Fc*(dSec-0.4*xPos)+Fs2*(dSec-cover);
         end

        clear auxMat auxVec auxVec1 auxVec2 xPos index as asR1 nRebar asMin...
            count M1 Med Ned x1 x2 xPos1 xPos2 space index1 tol Fs1 Fs2 Fc asBeam...
            iFrame bSec hSec dSec

        % column rebar   
        MrdRqCol(1,1)=1.3*MrdBeam(1)/2;
        MrdRqCol(2,1)=max([1.3*MrdBeam(1) 1.3*MrdBeam(2)/2]);
        MrdRqCol(3,1)=1.3*MrdBeam(2);

        calcAsCol=cell(1,size(frame,2));
        for iFrame=1:size(frame,2);
            count=0; auxMat=zeros(endColTag,1);
            for iElem=1:endColTag;
                 count=count+1; auxVec=zeros(nComb,1);
                 if iFrame==1;
                    if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                        hSec=secCol(1,1);bSec=secCol(1,2);
                    else
                        hSec=secCol(2,1);bSec=secCol(2,2);
                    end
                else
                    if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                        hSec=secCol(2,1);bSec=secCol(2,2);
                    else
                        hSec=secCol(3,1);bSec=secCol(3,2);
                    end
                 end
                 dSec=hSec-cover;
                 for iComb=1:nComb;
                     Med=max(abs(M{iComb,iFrame}(iElem,:)));
                     if iComb==8
                        if iFrame==1;
                            if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                                MrdRq=MrdRqCol(1);
                            else
                                MrdRq=MrdRqCol(2);
                            end
                        else
                            if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                                MrdRq=MrdRqCol(2);
                            else
                                MrdRq=MrdRqCol(3);
                            end
                        end
                       Med=max([Med MrdRq]);
                     end

                    Ned=N{iComb,iFrame}(iElem,1);
                    niu=Ned/(bSec*hSec*(fck*1000/1.5));
                    ecc= Med/Ned;
                    sigma1=Ned/(bSec*hSec)+Med/((bSec*hSec^3)/12)*(hSec/2);
                    sigma2=Ned/(bSec*hSec)-Med/((bSec*hSec^3)/12)*(hSec/2);
                    if sign(sigma1)==-1 && sign(sigma2)==-1; % simple compression
                        asTot=(abs(Ned)-(bSec*hSec*fcd))/400e3;
                    elseif sigma1<fctm*1000 && sign(sigma2)==-1; % column still in elastic range
                        asTot=(abs(Ned)-(bSec*hSec*fcd))/400e3;            
                    elseif abs(ecc)>hSec && sign(Ned)==-1;% design as regular bending with symmetrical rebar
                        tol=1;xPos=(cover-0.001); Mrd=0;
                        while tol>0.01&&xPos<hSec&&abs(Mrd)<Med;
                            xPos=xPos+0.001;
                            Fc=0.8*xPos*bSec*fck*1000/1.5;
                            aux1=xPos/(hSec-2*cover); aux2=((hSec-2*cover)-xPos)/(hSec-2*cover);
                            asCalc1=(abs(Ned)-Fc)/(200e6*(3.5e-3*(xPos-cover)/xPos)+...
                                2*aux1*200e6*(3.5e-3)/2-2*aux2*200e6*3.5e-3*((dSec-xPos)/2)/xPos-...
                                fy*1000/1.15);

                            Fs2=asCalc1*200e6*(3.5e-3*(xPos-cover)/xPos);
                            Fs3=2*asCalc1*aux1*200e6*(3.5e-3/2);
                            Fs4=2*asCalc1*aux2*200e6*(3.5e-3*((dSec-xPos)/2)/xPos);
                            Mrd=Fc*(dSec-0.4*xPos)+Fs2*(dSec-cover)+Fs3*(dSec-xPos)-Fs4*(dSec-xPos)/2;                    
                            tol=abs((Mrd-Med)/Med);
                        end
                        asTot=asCalc1*4;
                    end
                    auxVec(iComb,1)=max([0.10*abs(Ned)/(fy*1000/1.15) 0.002*bSec*hSec asTot]);            
                 end
                 auxMat(count,1)=max(auxVec);
            end
            calcAsCol{1,iFrame}=auxMat;
        end
        nPhi=[4;8;12;16];
        auxVec1=[];auxVec2=[];auxVec3=[];
        for iFrame=1:size(frame,2);
            count=0;
            for iElem=1:endColTag;
                count=count+1;
                if iFrame==1;
                    if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                        auxVec1=[auxVec1 calcAsCol{1,iFrame}(count,1)];
                    else
                        auxVec2=[auxVec2 calcAsCol{1,iFrame}(count,1)];
                    end            
                else
                     if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                         auxVec2=[auxVec2 calcAsCol{1,iFrame}(count,1)];
                     else
                        auxVec3=[auxVec3 calcAsCol{1,iFrame}(count,1)];
                     end
                end           
            end         
        end
        asCol(1,1)=max(auxVec1);asCol(2,1)=max(auxVec2);asCol(3,1)=max(auxVec3);
        rhoAsCol=asCol./(secCol(:,1).^2);

        lAsMat=zeros(size(phiRebar,1),size(nPhi,1));
        for iS=1:size(phiRebar,1);
            for iNPhi=1:size(nPhi,1);
                lAsMat(iS,iNPhi)=nPhi(iNPhi)*pi*phiRebar(iS)^2/4;
            end        
        end
        for iColSec=1:size(secCol,1);
            selectMat=lAsMat-(asCol(iColSec,1)*ones(size(phiRebar,1), size(nPhi,1)));
                auxVec=[];
            for iS=1:size(phiRebar,1);
                for iNPhi=1:size(nPhi,1);
                    if selectMat(iS,iNPhi)>0;
                        auxVec=[auxVec selectMat(iS,iNPhi)];
                    end
                end
            end
            minAux=min(auxVec);
           [index1,index2]=find(selectMat==minAux,1,'last');
           secCol(iColSec,3)=nPhi(index2);
           secCol(iColSec,4)=phiRebar(index1);
        end

        clear iFrame iElem auxMat auxVec Med Mrd auxVec1 auxVec2 auxVe3...
        aux1 aux2 Fc Fs Fs1 Fs2 Fs3 Fs4 xPos Ned sigma1 sigma2 bSec hSec dSec secDim...
        minAux index1 index2 lAsMat iS iNPhi
        %% design control
        miuCol(1,1)=(1.3*MrdBeam(1)/2)/(secCol(1,1)*secCol(1,2)^2*fck*1000/1.5);
        miuCol(2,1)=(max([1.3*MrdBeam(1) 1.3*MrdBeam(2)/2]))/(secCol(2,1)*secCol(2,2)^2*fck*1000/1.5);
        miuCol(3,1)=(1.3*MrdBeam(2))/(secCol(3,1)*secCol(3,2)^2*fck*1000/1.5);

        if dsgAccel>0;
            index=find(miuCol>0.25);
            controlDuct=isempty(index);
            if controlDuct==0;
                controlVal1=0;
                for iLine=index(1):index(end);
                    secDim=ceil(((MrdRqCol(iLine)/(0.25*fck*1000/1.5))^(1/3))/0.05)*0.05;
                    secCol(iLine,1)=secDim;
                    secCol(iLine,2)=secDim;
                end
            else
                controlVal1=1;
            end
        else
            controlVal1=1;
        end
        index=find(rhoAsCol>0.015);controlRhoL=isempty(index);
        if controlRhoL==0 && controlDuct==1;
            controlVal2=0;
             for iLine=index(1):index(end);
                 secCol(iLine,1)=secCol(iLine,1)+0.05;
                 secCol(iLine,2)=secCol(iLine,2)+0.05;
             end
        else
            controlVal2=1;
        end
        controlVal=min([controlVal1 controlVal2]);


        nIter=nIter+1;
        if nIter==10;
            error('Error: Script did not converge. Stable design could not be found.');
        end        
    end

    %% shear rebar
    sVec=[0.01 0.125 0.15 0.175 0.25]; phiVec=[0.006 0.008 0.010 0.012];
    aswMat=zeros(size(phiVec,2),size(sVec,2));

    for iPhi=1:size(phiVec,2);
        for iS=1:size(sVec,2);
            aswMat(iPhi,iS)=(2*pi*(phiVec(iPhi)^2)/4)/sVec(iS);
        end
    end

    % beam rebar
    for iFrame=1:size(frame,2)
        maxVBeam=max(max(VEnvAbsMax{1,iFrame}(endColTag+1:end,:)));
        hSec=secBeam(iFrame,2);bSec=secBeam(iFrame,1);
        dSec=hSec-cover; miu1=0.6*(1-fck/250);    
        zSec=0.9*dSec;
        cotTheta=roots([maxVBeam -1*(bSec*zSec*miu1*(fck*1000)/1.5) maxVBeam]);
        if max(cotTheta)>=2.5;
            cotTheta=2.5;
        elseif max(cotTheta)>=1
            cotTheta=max(cotTheta);
        else
            cotTheta=1;
        end
        aswOs=max([maxVBeam/(zSec*(fy*1000)/1.15*cotTheta) 0.08*sqrt(fck)/fy*bSec]);
        slMax=0.75*dSec*(1+(1/tan(pi/2)));
        selectMat=aswMat-aswOs;
        auxVec=[];
        for iPhi=1:size(phiVec,2);
            for iS=1:size(sVec,2);
                if selectMat(iPhi,iS)>0&&sVec(iS)<slMax;
                    auxVec=[auxVec selectMat(iPhi,iS)];
                end
            end
        end
        minAux=min(auxVec);
        [index1,index2]=find(selectMat==minAux);
        secBeam(iFrame,end-1)=phiVec(index1);
        secBeam(iFrame,end)=sVec(index2);    
    end

    %column rebar
    auxVec1=[];auxVec2=[];auxVec3=[];
    for iFrame=1:size(frame,2);
        count=0;
        for iElem=1:endColTag;
            count=count+1;
            if iFrame==1;
                if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                    auxVec1=[auxVec1 max(VEnvAbsMax{1,iFrame}(count,:))];
                else
                    auxVec2=[auxVec2 max(VEnvAbsMax{1,iFrame}(count,:))];
                end            
            else
                 if iElem<=nFloor||iElem>=(endColTag-nFloor+1);
                     auxVec2=[auxVec2 max(VEnvAbsMax{1,iFrame}(count,:))];
                 else
                     auxVec3=[auxVec3 max(VEnvAbsMax{1,iFrame}(count,:))];
                 end
            end           
        end
    end
    maxVCol(1,1)=max(auxVec1);maxVCol(2,1)=max(auxVec2);
    maxVCol(3,1)=max(auxVec3);
    for iLine=1:size(maxVCol,1);
        hSec=secCol(iLine,2);bSec=secCol(iLine,1);
        dSec=hSec-cover; miu1=0.6*(1-fck/250);    
        zSec=0.9*dSec;
        cotTheta=roots([maxVCol(iLine,1) -1*(bSec*zSec*miu1*(fck*1000)/1.5) maxVCol(iLine,1)]);
        if max(cotTheta)>=2.5;
            cotTheta=2.5;
        elseif max(cotTheta)>=1
            cotTheta=max(cotTheta);
        else
            cotTheta=1;
        end
        aswOs=max([maxVCol(iLine,1)/(zSec*(fy*1000)/1.15*cotTheta) 0.08*sqrt(fck)/fy*bSec]);
        slMax=0.6*min([20*secCol(iLine,4) min([secCol(iLine,1) secCol(iLine,2)]) 0.4]);
        phiMin=min([0.006 1/4*secCol(iLine,4)]);
        for iPhi=1:size(phiVec,2);
            for iS=1:size(sVec,2);
                if selectMat(iPhi,iS)>0&&sVec(iS)<slMax&&phiVec(iPhi)<phiMin;
                    auxVec=[auxVec selectMat(iPhi,iS)];
                end
            end
        end
        minAux=min(auxVec);
        [index1,index2]=find(selectMat==minAux);
        secCol(iLine,end-1)=phiVec(index1);
        secCol(iLine,end)=sVec(index2);        
    end
     clear iFrame iLine iS iPhi iElem auxMatAbsV auxVec auxVec1 auxVec2 auxVec3...
         minAux index1 index2
     %%
     % Compute confinement
     confCoefCol=zeros(size(secCol,1),1);
     confCoefBeam=zeros(size(secBeam,1),2); %[M- M+]
     % columns
    for iLine=1:size(secCol,1);
        bc=secCol(iLine,2)-(2*cover); dc=secCol(iLine,1)-(2*cover);
        phiT=secCol(iLine,end-1); sT=secCol(iLine,end);
        wiAvg=(2*bc+2*dc)/(secCol(iLine,3));
        rhoCC=(secCol(iLine,3)*pi*secCol(iLine,4)^2/4)/(bc*dc);
        ke=((1-(secCol(iLine,3)*(wiAvg^2)/6))*(1-sT/(2*bc))*(1-sT/(2*dc)))/...
            (1-rhoCC);
        rhoS=(pi*phiT^2/4)*(2*bc+2*dc)/(sT*bc*dc);
        fl2=1/2*ke*rhoS*fy;
        confCoefCol(iLine,1)=-1.254+2.254*sqrt(1+7.94*fl2/fck)-2*(fl2/fck);
    end

    % beams
     for iLine=1:size(secBeam,1);
          phiT=secBeam(iLine,end-1); sT=secBeam(iLine,end);
          bc=secBeam(iLine,2)-(2*cover); dc=secBeam(iLine,1)-(2*cover);
          rhoS=(pi*phiT^2/4)*(2*bc+2*dc)/(sT*bc*dc); 
          rhoCC1=(secBeam(iLine,3)*pi*secBeam(iLine,4)^2/4+...
              secBeam(iLine,5)*pi*secBeam(iLine,6)^2/4)/(bc*dc);
          rhoCC2=(secBeam(iLine,7)*pi*secBeam(iLine,8)^2/4+...
              secBeam(iLine,9)*pi*secBeam(iLine,10)^2/4)/(bc*dc);
          aux1=2*(bc^2)/6;aux2=2*(bc^2)/6;
          for iSide=1:2;
              nRebar1=secBeam(iLine,3+(iSide-1)*2);
              nRebar2=secBeam(iLine,7+(iSide-1)*2);
              if nRebar1==0;
                  aux1=aux1+dc^2/6;
              else
                  aux1=aux1+(dc/(nRebar2-1))^2/6;
              end
              if nRebar2==0;
                  aux2=aux2+dc^2/6;
              else
                  aux2=aux2+(dc/(nRebar2-1))^2/6;
              end              
          end
          ke1=((1-aux1)*(1-sT/(2*bc))*(1-sT/(2*dc)))/...
            (1-rhoCC1);
          ke2=((1-aux2)*(1-sT/(2*bc))*(1-sT/(2*dc)))/...
            (1-rhoCC2);
         fl21=1/2*ke1*rhoS*fy; fl22=1/2*ke2*rhoS*fy;
         confCoefBeam(iLine,1)=-1.254+2.254*sqrt(1+7.94*fl21/fck)-2*(fl21/fck);
         confCoefBeam(iLine,2)=-1.254+2.254*sqrt(1+7.94*fl22/fck)-2*(fl22/fck);
     end

     clear bc dc rhoS rhoCC rhoCC1 rhoCC2 aux1 aux2 nRebar1 nRebar2 ke ke1...
         ke2 iLine iSide fl2 fl21 fl22

    % ----------------------------------------------------------------------

     %% ---- write final 3D model ----
    fid=fopen('nonlinear3DModel.tcl','wt+');
    fprintf(fid,'# --Model Info--');fprintf(fid,'\n');
    fprintf(fid,'# nBay=%g; nFloor=%g; PGA[g]=%g \t',nBay,nFloor,dsgAccel/9.81);
    fprintf(fid,'\n');fprintf(fid,'# floor height[m]=%g; span[m]=%g \t',fHeight,spanL);
    fprintf(fid,'\n');fprintf(fid,'\n');
    fprintf(fid,'wipe all');
    fprintf(fid,'\n');fprintf(fid,'model BasicBuilder -ndm 3 -ndf 6');
    fprintf(fid,'\n');

    % -Nodes and Masses-
    fprintf(fid,'# --create nodes and masses--'); fprintf(fid,'\n');
    for iNode=1:nNodes3D;
        fprintf(fid,'node %g %g %g %g \t',nCoord3D(iNode,1),nCoord3D(iNode,2),...
            nCoord3D(iNode,3),nCoord3D(iNode,4));fprintf(fid,'\n');
        if nCoord3D(iNode,3)==0;
            fprintf(fid,'fix %g %g %g %g %g %g %g \t',nCoord3D(iNode,1),1,1,1,...
                1,1,1); fprintf(fid,'\n');
        end
        if mass3D(iNode,1)>0;
            fprintf(fid,'mass %g %g %g %g %g %g %g \t',nCoord3D(iNode,1),mass3D(iNode,1),...
            mass3D(iNode,1),mass3D(iNode,1),1e-9,1e-9,1e-9);fprintf(fid,'\n');    
        end
    end
    clear iNode

    % -Rigid diaphragm-
    fprintf(fid,'# --create rigid diaphragm--'); fprintf(fid,'\n');
    mNodeTag=9000;% master node tag
    for iFloor=1:nFloor;
        mNodeTag=mNodeTag+1;
        index=find(nCoord3D(:,3)==(fHeight*iFloor));
        % write masternode
        xMN=(nBay*spanL)/2; yMN=fHeight*iFloor; zMN=(nBay*spanL)/2;
        fprintf(fid,'node %g %g %g %g \t',mNodeTag,xMN,yMN,zMN); 
        fprintf(fid,'\n');
        fprintf(fid,'fix %g %g %g %g %g %g %g \t',mNodeTag,0,1,0,1,0,1);
        fprintf(fid,'\n');
        % write rigid diaphragm
        sNodeTag=nCoord3D(index,1); % slave nodes' tag
        fprintf(fid,'rigidDiaphragm %g %g \t',2,mNodeTag);
        for iLine=1:length(sNodeTag);
            fprintf(fid,'%g \t',sNodeTag(iLine));
        end
        fprintf(fid,'\n');       
    end
    clear iFloor index xMN yMN xMN iLine

    % -Materials-
    fprintf(fid,'# --create uniaxialMaterial--'); fprintf(fid,'\n');
    steelTag=100;uConcrTag=200;cConcrTag=300;
    % steel
    fprintf(fid,'uniaxialMaterial Steel02 %g %g %g %g %g %g %g \t',steelTag,...
        fy*1000,200e6,0.005,18,0.925,0.15); fprintf(fid,'\n');
    % unconfined concrete   
    fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
        uConcrTag,-1*fcm*1000,-1*epsc0,-1*epscU,5000*sqrt(fcm)*1000); 
    fprintf(fid,'\n');
    % confined concrete
    % ---columns---
    for iLine=1:size(confCoefCol,1);
        cConcrTag=cConcrTag+1;
        fcc=-1*fcm*confCoefCol(iLine,1)*1000;
        epscc0=-1*epsc0*(1+5*(confCoefCol(iLine,1)-1));
        epsccU=-1*epscU*(1+5*(confCoefCol(iLine,1)-1));
        fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
            cConcrTag,fcc,epscc0,epsccU,5000*sqrt(fcm)*1000);
        fprintf(fid,'\n');    
    end
    % ---beams---
    for iLine=1:size(confCoefBeam,1);
        for iCol=1:size(confCoefBeam,2);
            cConcrTag=cConcrTag+1;
            fcc=-1*fcm*confCoefBeam(iLine,iCol)*1000;
            epscc0=-1*epsc0*(1+5*(confCoefBeam(iLine,iCol)-1));
            epsccU=-1*epscU*(1+5*(confCoefBeam(iLine,iCol)-1));
            fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
                cConcrTag,fcc,epscc0,epsccU,5000*sqrt(fcm)*1000);
            fprintf(fid,'\n');
        end
    end
    clear iLine iCol fcc epsccO epsccU

    % -fiber sections- 
    fprintf(fid,'# --create fiber sections--'); fprintf(fid,'\n');
    fSecTag=0; cConcrTag=300;
    JxVec=[];
    % columns
    secC=[];
    for iLine=1:size(secCol,1);
        fSecTag=fSecTag+1; cConcrTag=cConcrTag+1;
        secC=[secC fSecTag];
        bw=secCol(iLine,1); hw=secCol(iLine,2);
        y1=hw/2;
        z1=bw/2;
        nRebar=secCol(iLine,3);phiRebar=secCol(iLine,4);
        if y1>=z1;
            a=y1;b=z1;
        else
            a=z1;b=y1;
        end
        JxVec=[JxVec a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))]; 
        % start section fiber command
        fprintf(fid,'section Fiber %g { \t',fSecTag);
        fprintf(fid,'\n');
        % core concrete fibers
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',cConcrTag,10,...
            10,(cover-y1),(cover-z1),(y1-cover),(z1-cover));
        fprintf(fid,'\n');
        % unconfined concrete fibers
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,z1-cover,y1,z1);
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,-z1,y1,(cover-z1));
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,(cover-z1),(cover-y1),(z1-cover));
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,(y1-cover),(cover-z1),y1,(z1-cover));
        fprintf(fid,'\n');
        % steel fibers
        % corner rebar
        fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
            2,pi*phiRebar^2/4,(y1-cover),(z1-cover),(y1-cover),...
            (cover-z1));
        fprintf(fid,'\n');
        fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
            2,pi*phiRebar^2/4,(cover-y1),(z1-cover),(cover-y1),...
            (cover-z1));
        fprintf(fid,'\n');
        if nRebar==8;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover),0,(cover-y1),0);
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,0,(z1-cover),0,(cover-z1));
            fprintf(fid,'\n');
        elseif nRebar==12
            space1=(bw-2*cover)/3;space2=(hw-2*cover)/3;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover-space2),(z1-cover),(y1-cover-2*space2),(z1-cover));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover-space2),(cover-z1),(y1-cover-2*space2),(cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover),(z1-cover-space1),(y1-cover),(z1-cover-2*space1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(cover-y1),(z1-cover-space1),(cover-y1),(z1-cover-2*space1));
            fprintf(fid,'\n');
        elseif nRebar==16;
            space1=(bw-2*cover)/4;space2=(hw-2*cover)/4;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover-space2),(z1-cover),(y1-cover-2*space2),(z1-cover));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover-space2),(cover-z1),(y1-cover-2*space2),(cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover),(z1-cover-space1),(y1-cover),(z1-cover-2*space1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(cover-y1),(z1-cover-space1),(cover-y1),(z1-cover-2*space1));
            fprintf(fid,'\n');        
        end    
        % close section fiber command
        fprintf(fid,'} \t');
        fprintf(fid,'\n');   
    end
    % beams
    secB=[];
    for iLine=1:size(confCoefBeam,1);
        for iCol=1:size(confCoefBeam,2);
            fSecTag=fSecTag+1; cConcrTag=cConcrTag+1;
            secB=[secB fSecTag];
            bw=secBeam(iLine,1); hw=secBeam(iLine,2);
            y1=hw/2;
            z1=bw/2;
            nRebar=secCol(iLine,3);phiRebar=secCol(iLine,4);
            if y1>=z1;
                a=y1;b=z1;
            else
                a=z1;b=y1;
            end
            JxVec=[JxVec a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))];
            % start section fiber command
            fprintf(fid,'section Fiber %g { \t',fSecTag);
            fprintf(fid,'\n');
            % core concrete fibers
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',cConcrTag,10,...
                10,(cover-y1),(cover-z1),(y1-cover),(z1-cover));
            fprintf(fid,'\n');
            % unconfined concrete fibers
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,z1-cover,y1,z1);
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,-z1,y1,(cover-z1));
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,(cover-z1),(cover-y1),(z1-cover));
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,(y1-cover),(cover-z1),y1,(z1-cover));
            fprintf(fid,'\n');
            % steel fibers
            nRebarUp=secBeam(iLine,3+4*(iCol-1));
            phiRebarUp=secBeam(iLine,4+4*(iCol-1));
            nRebarDwn=secBeam(iLine,5+4*(iCol-1));
            phiRebarDwn=secBeam(iLine,6+4*(iCol-1));
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                nRebarUp,pi*phiRebarUp^2/4,(y1-cover),(z1-cover),(y1-cover),...
                    (cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                nRebarDwn,pi*phiRebarDwn^2/4,(cover-y1),(z1-cover),(cover-y1),...
                    (cover-z1));
            fprintf(fid,'\n');        
            % close section fiber command
            fprintf(fid,'} \t');
            fprintf(fid,'\n');   

        end
    end
    nSec=fSecTag;
    clear iLine iCol bw hw y1 z1 a b

    % -aggregate torsion-
    fprintf(fid,'# --aggregate torsion--'); fprintf(fid,'\n');
    aSecTag=900; % aggregated section tag;
    tMatTag=900;
    for iSec=1:nSec
        aSecTag=aSecTag+1; tMatTag=tMatTag+1;
        fprintf(fid,'uniaxialMaterial Elastic %g %g \t',tMatTag,Gc*1000*JxVec(iSec));
        fprintf(fid,'\n');
        fprintf(fid,'section Aggregator %g %g T -section %g \t',...
                 aSecTag,tMatTag,iSec);
        fprintf(fid,'\n');    
    end
    aSecTag=900;
    % -geomTransf- 
    fprintf(fid,'# --create geomTransf--'); fprintf(fid,'\n');
    %columns
    fprintf(fid,'geomTransf Corotational %g %g %g %g\t',1,0,0,1); fprintf(fid,'\n');
    %beams
    fprintf(fid,'geomTransf Linear %g %g %g %g\t',2,0,0,1); fprintf(fid,'\n');
    %girders
    fprintf(fid,'geomTransf Linear %g %g %g %g\t',3,-1,0,0); fprintf(fid,'\n');

    % -elements-
    fprintf(fid,'# --create elements--'); fprintf(fid,'\n');
    % columns
    eleTag=0;
    for iBayZ=1:nBay+1;
        for iBayX=1:nBay+1;
            for iFloor=1:nFloor;
                eleTag=eleTag+1;
                xCoordI=(iBayX-1)*spanL;yCoordI=(iFloor-1)*fHeight;
                zCoordI=(iBayZ-1)*spanL;
                xCoordJ=(iBayX-1)*spanL;yCoordJ=(iFloor)*fHeight;
                zCoordJ=(iBayZ-1)*spanL;
                auxMat1=nCoord3D(nCoord3D(:,2)==xCoordI,:);
                auxMat1=auxMat1(auxMat1(:,3)==yCoordI,:);
                nodeI=auxMat1(auxMat1(:,4)==zCoordI);
                auxMat2=nCoord3D(nCoord3D(:,2)==xCoordJ,:);
                auxMat2=auxMat2(auxMat2(:,3)==yCoordJ,:);
                nodeJ=auxMat2(auxMat2(:,4)==zCoordJ);
                if iBayZ==1||iBayZ==nBay+1;
                    if iBayX==1||iBayX==nBay+1;
                        fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                            eleTag,nodeI,nodeJ,5,aSecTag+secC(1),1,50,1e-10);
                        fprintf(fid,'\n');
                    else
                        fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                            eleTag,nodeI,nodeJ,5,aSecTag+secC(2),1,50,1e-10);
                        fprintf(fid,'\n');
                    end
                else
                    if iBayX==1||iBayX==nBay+1;
                        fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                            eleTag,nodeI,nodeJ,5,aSecTag+secC(2),1,50,1e-10);
                        fprintf(fid,'\n');
                    else
                        fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                            eleTag,nodeI,nodeJ,5,aSecTag+secC(3),1,50,1e-10);
                        fprintf(fid,'\n');
                    end
                end            
            end
        end
    end
    endColTag3D=eleTag;
    % beams
    for iBayZ=1:nBay+1;
        for iBayX=1:nBay;
            for iFloor=1:nFloor;
                eleTag=eleTag+1;
                xCoordI=(iBayX-1)*spanL;yCoordI=(iFloor)*fHeight;
                zCoordI=(iBayZ-1)*spanL;
                xCoordJ=(iBayX)*spanL;yCoordJ=(iFloor)*fHeight;
                zCoordJ=(iBayZ-1)*spanL;
                auxMat1=nCoord3D(nCoord3D(:,2)==xCoordI,:);
                auxMat1=auxMat1(auxMat1(:,3)==yCoordI,:);
                nodeI=auxMat1(auxMat1(:,4)==zCoordI);
                auxMat2=nCoord3D(nCoord3D(:,2)==xCoordJ,:);
                auxMat2=auxMat2(auxMat2(:,3)==yCoordJ,:);
                nodeJ=auxMat2(auxMat2(:,4)==zCoordJ);
                if iBayZ==1||iBayZ==nBay+1;
                    fprintf(fid,'element forceBeamColumn %g %g %g %g -sections %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,aSecTag+secB(1),aSecTag+secB(1),...
                    aSecTag+secB(2),aSecTag+secB(1),aSecTag+secB(1),2,50,1e-10);
                    fprintf(fid,'\n');
                else
                    fprintf(fid,'element forceBeamColumn %g %g %g %g -sections %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,aSecTag+secB(3),aSecTag+secB(3),...
                    aSecTag+secB(4),aSecTag+secB(3),aSecTag+secB(3),2,50,1e-10);
                    fprintf(fid,'\n');
                end

            end
        end
    end
    %girders
    for iBayX=1:nBay+1;
        for iBayZ=1:nBay;
            for iFloor=1:nFloor;
                eleTag=eleTag+1;
                xCoordI=(iBayX-1)*spanL;yCoordI=(iFloor)*fHeight;
                zCoordI=(iBayZ-1)*spanL;
                xCoordJ=(iBayX-1)*spanL;yCoordJ=(iFloor)*fHeight;
                zCoordJ=(iBayZ)*spanL;
                auxMat1=nCoord3D(nCoord3D(:,2)==xCoordI,:);
                auxMat1=auxMat1(auxMat1(:,3)==yCoordI,:);
                nodeI=auxMat1(auxMat1(:,4)==zCoordI);
                auxMat2=nCoord3D(nCoord3D(:,2)==xCoordJ,:);
                auxMat2=auxMat2(auxMat2(:,3)==yCoordJ,:);
                nodeJ=auxMat2(auxMat2(:,4)==zCoordJ);
                if iBayX==1||iBayX==nBay+1;
                    fprintf(fid,'element forceBeamColumn %g %g %g %g -sections %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,aSecTag+secB(1),aSecTag+secB(1),...
                    aSecTag+secB(2),aSecTag+secB(1),aSecTag+secB(1),3,50,1e-10);
                    fprintf(fid,'\n');
                else
                    fprintf(fid,'element forceBeamColumn %g %g %g %g -sections %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,aSecTag+secB(3),aSecTag+secB(3),...
                   aSecTag+secB(4),aSecTag+secB(3),aSecTag+secB(3),3,50,1e-10);
                    fprintf(fid,'\n');
                end          
            end
        end
    end
    nElem3D=eleTag;
    clear nodeI nodeJ auxMat1 auxMat2 iFloor iBayX iBayZ xCoordI xCoordJ yCoordI yCoordJ...
        zCoordI zCoordJ
    
    % eigen analysis
    fprintf(fid,'# --eigen analysis--'); fprintf(fid,'\n');
    fprintf(fid,'set pi %g',pi);fprintf(fid,'\n');
    fprintf(fid,'set lambda [eigen -fullGenLapack 1]');fprintf(fid,'\n');
    fprintf(fid,'set T [expr 2.0*$pi/pow([lindex $lambda 0],0.5)]');fprintf(fid,'\n');
    fprintf(fid,'set outfile [open "periodEl.txt" w]');fprintf(fid,'\n');
    fprintf(fid,'puts $outfile $T');fprintf(fid,'\n');
    fprintf(fid,'recorder Node -file phi1.txt -nodeRange 9001 %g -dof 1 "eigen 1"',mNodeTag);
    fprintf(fid,'\n');
    % -gravity loads
    fprintf(fid,'# --apply gravity loads--'); fprintf(fid,'\n');
    fprintf(fid,'pattern Plain 1 Linear {');fprintf(fid,'\n');
    eleTag=endColTag3D; bLoad3D=zeros(nElem3D,1);
    for iBayZ=1:nBay+1;
        for iBayX=1:nBay;
            for iFloor=1:nFloor; 
                eleTag=eleTag+1;
                if iBayZ==1||iBayZ==nBay+1;
                    iflAreaB=(spanL*(spanL/2))/2;
                else
                     iflAreaB=((spanL*(spanL/2))/2)*2;
                end
                if iFloor==nFloor
                    live=qkTop;
                else 
                    live=qk;
                end
                     bLoad3D(eleTag,1)=(-1*(gk+psi2(2)*live)*iflAreaB)/spanL;
                     fprintf(fid,'eleLoad -ele %g -type -beamUniform %g %g \t',...
                             eleTag,bLoad3D(eleTag,1),0); fprintf(fid,'\n');         
            end
        end
    end
    for iBayX=1:nBay+1;
        for iBayZ=1:nBay;
            for iFloor=1:nFloor; 
                eleTag=eleTag+1;
                if iBayX==1||iBayX==nBay+1;
                    iflAreaB=(spanL*(spanL/2))/2;
                else
                     iflAreaB=((spanL*(spanL/2))/2)*2;
                end
                if iFloor==nFloor
                    live=qkTop;
                else 
                    live=qk;
                end
                     bLoad3D(eleTag,1)=(-1*(gk+psi2(2)*live)*iflAreaB)/spanL;
                     fprintf(fid,'eleLoad -ele %g -type -beamUniform %g %g \t',...
                             eleTag,bLoad3D(eleTag,1),0); fprintf(fid,'\n');         
            end
        end
    end
    fprintf(fid,'} \t');fprintf(fid,'\n');

    fprintf(fid,'# --solve gravity loads--'); fprintf(fid,'\n');
    fprintf(fid,'set nSteps 21');fprintf(fid,'\n');
    fprintf(fid,'for {set iStep 1} {$iStep<$nSteps} {incr iStep 1} {');
    fprintf(fid,'\n');fprintf(fid,'system BandGeneral');
    fprintf(fid,'\n'); fprintf(fid,'constraints Transformation');
    fprintf(fid,'\n');fprintf(fid,'numberer RCM');fprintf(fid,'\n');
    fprintf(fid,'test NormDispIncr 1.0e-5  10 3');fprintf(fid,'\n');
    fprintf(fid,'algorithm Newton');fprintf(fid,'\n');fprintf(fid,'integrator LoadControl 0.05');
    fprintf(fid,'\n');fprintf(fid,'analysis Static');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]');fprintf(fid,'\n');
    fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
    fprintf(fid,'\n'); fprintf(fid,'algorithm Newton -initial');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]}');fprintf(fid,'\n');
    fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
    fprintf(fid,'\n'); fprintf(fid,'algorithm ModifiedNewton -initial');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]}}');fprintf(fid,'\n');
    fprintf(fid,'loadConst -time 0.0');fprintf(fid,'\n');

    clear iBayZ iBayx iFloor iflAreaB live secC secB

    % close model file
    fclose(fid);

    %% ---- write final 2D model (interior frame) ----
    fid=fopen('nonlinear2DModel.tcl','wt+');
    fprintf(fid,'# --Model Info--');fprintf(fid,'\n');
    fprintf(fid,'# nBay=%g; nFloor=%g; PGA[g]=%g \t',nBay,nFloor,dsgAccel/9.81);
    fprintf(fid,'\n');fprintf(fid,'# floor height[m]=%g; span[m]=%g \t',fHeight,spanL);
    fprintf(fid,'\n');fprintf(fid,'\n');
    fprintf(fid,'wipe all');
    fprintf(fid,'\n');fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3');
    fprintf(fid,'\n');

    % -Nodes and Masses-
    fprintf(fid,'# --create nodes and masses--'); fprintf(fid,'\n');
    for iNode=1:nNodes2D;
        fprintf(fid,'node %g %g %g \t',nCoord2D(iNode,1),nCoord2D(iNode,2),...
            nCoord2D(iNode,3));fprintf(fid,'\n');
        if nCoord2D(iNode,3)==0;
            fprintf(fid,'fix %g %g %g %g \t',nCoord2D(iNode,1),1,1,1);
            fprintf(fid,'\n');
        end
        if mass2D(iNode,2)>0;
            fprintf(fid,'mass %g %g %g %g \t',nCoord2D(iNode,1),mass2D(iNode,2),...
            mass2D(iNode,2),1e-9);fprintf(fid,'\n');    
        end
    end
    clear iNode

    % -Materials-
    fprintf(fid,'# --create uniaxialMaterial--'); fprintf(fid,'\n');
    steelTag=100;uConcrTag=200;cConcrTag=300;
    % steel
    fprintf(fid,'uniaxialMaterial Steel02 %g %g %g %g %g %g %g \t',steelTag,...
        fy*1000,200e6,0.005,18,0.925,0.15); fprintf(fid,'\n');
    % unconfined concrete   
    fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
        uConcrTag,-1*fcm*1000,-1*epsc0,-1*epscU,5000*sqrt(fcm)*1000); 
    fprintf(fid,'\n');
    % confined concrete
    % ---columns---
    for iLine=2:size(confCoefCol,1);
        cConcrTag=cConcrTag+1;
        fcc=-1*fcm*confCoefCol(iLine,1)*1000;
        epscc0=-1*epsc0*(1+5*(confCoefCol(iLine,1)-1));
        epsccU=-1*epscU*(1+5*(confCoefCol(iLine,1)-1));
        fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
            cConcrTag,fcc,epscc0,epsccU,5000*sqrt(fcm)*1000);
        fprintf(fid,'\n');    
    end
    % ---beams---
    for iLine=2:size(confCoefBeam,1);
        for iCol=1:size(confCoefBeam,2);
            cConcrTag=cConcrTag+1;
            fcc=-1*fcm*confCoefBeam(iLine,iCol)*1000;
            epscc0=-1*epsc0*(1+5*(confCoefBeam(iLine,iCol)-1));
            epsccU=-1*epscU*(1+5*(confCoefBeam(iLine,iCol)-1));
            fprintf(fid,'uniaxialMaterial Concrete04 %g %g %g %g %g \t',...
                cConcrTag,fcc,epscc0,epsccU,5000*sqrt(fcm)*1000);
            fprintf(fid,'\n');
        end
    end
    clear iLine iCol fcc epsccO epsccU

    % -fiber sections- 
    fprintf(fid,'# --create fiber sections--'); fprintf(fid,'\n');
    fSecTag=0; cConcrTag=300;secC=[];
    JxVec=[];
    % columns
    for iLine=2:size(secCol,1);
        fSecTag=fSecTag+1; cConcrTag=cConcrTag+1;
        secC=[secC fSecTag];
        bw=secCol(iLine,1); hw=secCol(iLine,2);
        y1=hw/2;
        z1=bw/2;
        nRebar=secCol(iLine,3);phiRebar=secCol(iLine,4);
        if y1>z1;
            a=y1;b=z1;
        else
            a=z1;b=y1;
        end
        JxVec=[JxVec a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))]; 
        % start section fiber command
        fprintf(fid,'section Fiber %g { \t',fSecTag);
        fprintf(fid,'\n');
        % core concrete fibers
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',cConcrTag,10,...
            10,(cover-y1),(cover-z1),(y1-cover),(z1-cover));
        fprintf(fid,'\n');
        % unconfined concrete fibers
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,z1-cover,y1,z1);
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,-z1,y1,(cover-z1));
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,-y1,(cover-z1),(cover-y1),(z1-cover));
        fprintf(fid,'\n');    
        fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
            8,(y1-cover),(cover-z1),y1,(z1-cover));
        fprintf(fid,'\n');
        % steel fibers
        % corner rebar
        fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
            2,pi*phiRebar^2/4,(y1-cover),(z1-cover),(y1-cover),...
            (cover-z1));
        fprintf(fid,'\n');
        fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
            2,pi*phiRebar^2/4,(cover-y1),(z1-cover),(cover-y1),...
            (cover-z1));
        fprintf(fid,'\n');
        if nRebar==8;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover),0,(cover-y1),0);
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,0,(z1-cover),0,(cover-z1));
            fprintf(fid,'\n');
        elseif nRebar==12
            space1=(bw-2*cover)/3;space2=(hw-2*cover)/3;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover-space2),(z1-cover),(y1-cover-2*space2),(z1-cover));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover-space2),(cover-z1),(y1-cover-2*space2),(cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(y1-cover),(z1-cover-space1),(y1-cover),(z1-cover-2*space1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                2,pi*phiRebar^2/4,(cover-y1),(z1-cover-space1),(cover-y1),(z1-cover-2*space1));
            fprintf(fid,'\n');
        elseif nRebar==16;
            space1=(bw-2*cover)/4;space2=(hw-2*cover)/4;
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover-space2),(z1-cover),(y1-cover-2*space2),(z1-cover));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover-space2),(cover-z1),(y1-cover-2*space2),(cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(y1-cover),(z1-cover-space1),(y1-cover),(z1-cover-2*space1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                3,pi*phiRebar^2/4,(cover-y1),(z1-cover-space1),(cover-y1),(z1-cover-2*space1));
            fprintf(fid,'\n');        
        end        
        % close section fiber command
        fprintf(fid,'} \t');
        fprintf(fid,'\n');

    end
    % beams
    secB=[];
    for iLine=2:size(confCoefBeam,1);
        for iCol=1:size(confCoefBeam,2);
            fSecTag=fSecTag+1; cConcrTag=cConcrTag+1;
            secB=[secB fSecTag];
            bw=secBeam(iLine,1); hw=secBeam(iLine,2);
            y1=hw/2;
            z1=bw/2;
            if y1>z1;
                a=y1;b=z1;
            else
                a=z1;b=y1;
            end
            JxVec=[JxVec a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))];
            % start section fiber command
            fprintf(fid,'section Fiber %g { \t',fSecTag);
            fprintf(fid,'\n');
            % core concrete fibers
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',cConcrTag,10,...
                10,(cover-y1),(cover-z1),(y1-cover),(z1-cover));
            fprintf(fid,'\n');
            % unconfined concrete fibers
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,z1-cover,y1,z1);
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,-z1,y1,(cover-z1));
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,-y1,(cover-z1),(cover-y1),(z1-cover));
            fprintf(fid,'\n');    
            fprintf(fid,'patch rect %g %g %g %g %g %g %g \t',uConcrTag,8,...
                8,(y1-cover),(cover-z1),y1,(z1-cover));
            fprintf(fid,'\n');
            % steel fibers
            nRebarUp=secBeam(iLine,3+4*(iCol-1));
            phiRebarUp=secBeam(iLine,4+4*(iCol-1));
            nRebarDwn=secBeam(iLine,5+4*(iCol-1));
            phiRebarDwn=secBeam(iLine,6+4*(iCol-1));
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                nRebarUp,pi*phiRebarUp^2/4,(y1-cover),(z1-cover),(y1-cover),...
                    (cover-z1));
            fprintf(fid,'\n');
            fprintf(fid,'layer straight %g %g %g %g %g %g %g \t',steelTag,...
                nRebarDwn,pi*phiRebarDwn^2/4,(cover-y1),(z1-cover),(cover-y1),...
                    (cover-z1));
            fprintf(fid,'\n');        
            % close section fiber command
            fprintf(fid,'} \t');
            fprintf(fid,'\n');        
        end
    end
    nSec=fSecTag;
    clear iLine iCol bw hw y1 z1 a b

    % -geomTransf- 
    fprintf(fid,'# --create geomTransf--'); fprintf(fid,'\n');
    fprintf(fid,'geomTransf Corotational %g \t',1);
    fprintf(fid,'\n');
    fprintf(fid,'geomTransf Linear %g \t',2);
    fprintf(fid,'\n');

    % -elements-
    fprintf(fid,'# --create elements--'); fprintf(fid,'\n');
    eleTag=0;
    % columns
    for iBay=1:nBay+1;
        for iFloor=1:nFloor;
            eleTag=eleTag+1;
            xCoordI=(iBay-1)*spanL;yCoordI=(iFloor-1)*fHeight;
            xCoordJ=(iBay-1)*spanL;yCoordJ=(iFloor)*fHeight;
            auxMat1=nCoord2D(nCoord2D(:,2)==xCoordI,:);
            auxMat2=nCoord2D(nCoord2D(:,2)==xCoordJ,:);
            nodeI=auxMat1((find(auxMat1(:,3)==yCoordI)),1);
            nodeJ=auxMat2((find(auxMat2(:,3)==yCoordJ)),1);
            if iBay==1||iBay==nBay+1
                fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,secC(1),1,50,1e-8);
                fprintf(fid,'\n');
            else
                fprintf(fid,'element forceBeamColumn %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,secC(2),1,50,1e-8);
                fprintf(fid,'\n');
            end          
        end
    end
    % beams
    for iFloor=1:nFloor;
        for iBay=1:nBay;
            eleTag=eleTag+1;
            xCoordI=(iBay-1)*spanL;yCoordI=(iFloor)*fHeight;
            xCoordJ=(iBay)*spanL;yCoordJ=(iFloor)*fHeight;
            auxMat1=nCoord2D(nCoord2D(:,2)==xCoordI,:);
            auxMat2=nCoord2D(nCoord2D(:,2)==xCoordJ,:);
            nodeI=auxMat1((find(auxMat1(:,3)==yCoordI)),1);
            nodeJ=auxMat2((find(auxMat2(:,3)==yCoordJ)),1);
            fprintf(fid,'element forceBeamColumn %g %g %g %g -sections %g %g %g %g %g %g -iter %g %g \t',...
                    eleTag,nodeI,nodeJ,5,secB(1),secB(1),secB(2),secB(1),secB(1),2,50,1e-8);
                fprintf(fid,'\n');
        end
    end
    clear nodeI nodeJ auxMat1 auxMat2 iFloor iBay xCoordI xCoordj yCoordI yCoordJ
    
    % eigen analysis
    fprintf(fid,'# --eigen analysis--'); fprintf(fid,'\n');
    fprintf(fid,'set pi %g',pi);fprintf(fid,'\n');
    fprintf(fid,'set lambda [eigen -fullGenLapack 1]');fprintf(fid,'\n');
    fprintf(fid,'set T [expr 2.0*$pi/pow([lindex $lambda 0],0.5)]');fprintf(fid,'\n');
    fprintf(fid,'set outfile [open "periodEl.txt" w]');fprintf(fid,'\n');
    fprintf(fid,'puts $outfile $T');fprintf(fid,'\n');
    ctrlNode1=min(nCoord2D(nCoord2D(:,3)==fHeight,1));
    ctrlNode2=min(nCoord2D(nCoord2D(:,3)==(fHeight*nFloor),1));
    fprintf(fid,'recorder Node -file phi1.txt -nodeRange %g %g -dof 1 "eigen1" \t',...
                ctrlNode1,ctrlNode2);fprintf(fid,'\n');

    % -gravity loads-
    fprintf(fid,'# --apply gravity loads--'); fprintf(fid,'\n');
    fprintf(fid,'pattern Plain 1 Linear {');fprintf(fid,'\n');
    eleTag=endColTag;
    bLoad=zeros(nElem,1);
    for iFloor=1:nFloor;
         for iBay=1:nBay;
             eleTag=eleTag+1;
             iflAreaB=((spanL*(spanL/2))/2)*2;
             if iFloor==nFloor
                 live=qkTop;
             else 
                 live=qk;
             end
             bLoad(eleTag,1)=(-1*(gk+psi2(2)*live)*iflAreaB)/spanL;
             fprintf(fid,'eleLoad -ele %g -type -beamUniform %g \t',...
                     eleTag,bLoad(eleTag,1)); fprintf(fid,'\n');
         end
    end
    fprintf(fid,'} \t'); fprintf(fid,'\n');

    fprintf(fid,'# --solve gravity loads--'); fprintf(fid,'\n');
    fprintf(fid,'set nSteps 21');fprintf(fid,'\n');
    fprintf(fid,'for {set iStep 1} {$iStep<$nSteps} {incr iStep 1} {');
    fprintf(fid,'\n');fprintf(fid,'system BandGeneral');
    fprintf(fid,'\n'); fprintf(fid,'constraints Transformation');
    fprintf(fid,'\n');fprintf(fid,'numberer RCM');fprintf(fid,'\n');
    fprintf(fid,'test NormDispIncr 1.0e-5  10 3');fprintf(fid,'\n');
    fprintf(fid,'algorithm Newton');fprintf(fid,'\n');fprintf(fid,'integrator LoadControl 0.05');
    fprintf(fid,'\n');fprintf(fid,'analysis Static');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]');fprintf(fid,'\n');
    fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
    fprintf(fid,'\n'); fprintf(fid,'algorithm Newton -initial');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]}');fprintf(fid,'\n');
    fprintf(fid,'if {$ok != 0} {');fprintf(fid,'\n'); fprintf(fid,'test NormDispIncr 1.0e-5 500 0');
    fprintf(fid,'\n'); fprintf(fid,'algorithm ModifiedNewton -initial');fprintf(fid,'\n');
    fprintf(fid,'set ok [analyze 1]}}');fprintf(fid,'\n');
    fprintf(fid,'loadConst -time 0.0');fprintf(fid,'\n');

    % close model file
    fclose(fid);
    %%
    rmdir('OutFrame2DExt','s');rmdir('OutFrame2DInt','s');
    delete('frame2DExt.tcl','frame2DInt.tcl');
    disp('Nonlinear models successfully created.')
    outFile=['modelData','nBay',num2str(nBay),...
        'nFloor',num2str(nFloor),'PGA',num2str(dsgAccel/9.81),'.mat'];
    save(char(outFile));

end