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