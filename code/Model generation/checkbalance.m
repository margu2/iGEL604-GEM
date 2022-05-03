
modelTemp=model1;
biomassreactions={'EX_','import','Teichuronic_acid','CellWall','aa_to_prot','CellMembrane','RNA','DNA','Glycogen','Cofactor_Pool','Biomass','EXC_OUT_Biomass'};
[modelTemp, metFormulae, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(modelTemp,'nameCMs',1);
[EleBal2, element2, metEle2] = checkEleBalance(modelTemp); %Calculates elemental balance for every reaction in the model

unbalrxnsC=[];
for i=1:length(modelTemp.rxns); 
    if EleBal2(find(strcmp(element2, 'C')),i)~=0 && ~contains(modelTemp.rxns(i), biomassreactions)
        unbalrxnsC=[unbalrxnsC; modelTemp.rxns(i), EleBal2(find(strcmp(element2, 'C')),i)];
        %Save reactions for which carbon is unbalanced and subsequently
        %delete them. (Excluding exchange and biomass reactions)
    end
end

unbalrxnsO=[];
for i=1:length(modelTemp.rxns); 
    if EleBal2(find(strcmp(element2, 'O')),i)~=0 && ~contains(modelTemp.rxns(i), biomassreactions)
        unbalrxnsO=[unbalrxnsO; modelTemp.rxns(i), EleBal2(find(strcmp(element2, 'O')),i)];
        %Save reactions for which oxygen is unbalanced and subsequently
        %delete them. (Excluding exchange and biomass reactions)
    end
end
modelTemp=removeRxns(modelTemp, unbalrxnsC(:,1), 'metFlag', false);


%Reperform elemental balance check since number of reactions has changed.  
[EleBal3, element3, metEle3] = checkEleBalance(modelTemp);
unbalrxnsH=[];
% Due to the incorrect metabolite formulas of KEGG (Metabolites are always fully
% protonated), the reactions often exclude protons on one of the sides of
% the reactions. Since the metabolites formulas are changed to MetaCyc
% formulas, these incorrect reactions need to be fixed. This script finds
% these reactions and automatically balance them. 
for i=1:length(modelTemp.rxns); 
    if EleBal3(find(strcmp(element2, 'H')),i)~=0 && ~contains(modelTemp.rxns(i), biomassreactions)
        unbalrxnsH=[unbalrxnsH; modelTemp.rxns(i), EleBal3(find(strcmp(element2, 'H')),i)];
        [metList, stoichiometries] = findMetsFromRxns(modelTemp, modelTemp.rxns(i));
        ADDRXNbal.mets=[metList{1}; 'C00080']; %Finds reactions unbalanced in hydrogen and adds a proton to the reaction side that needs it
        ADDRXNbal.stoichCoeffs=[stoichiometries{1}', -EleBal3(1,i)]; 
        modelTemp=changeRxns(modelTemp, modelTemp.rxns(i), ADDRXNbal, 1, 's', false);
    end
end
 
clear EleBal element metEle EleBal2 element2 metEle2 ADDRXNbal stoichiometries i
model1=modelTemp;