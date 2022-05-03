rxns={};
A=[];
modelTemp=model1;
for i=1:length(missingenesmodel.genes)
    if strcmp(modelTemp.genes, missingenesmodel.genes(i))~=1
        A=[A; missingenesmodel.genes(i)];
    end
end

for i=1:length(A)
    rxns(i)={missingenesmodel.rxns(find(contains(missingenesmodel.grRules, A{i})))};

end
rxns=rxns';

for i=1:length(rxns)
    if ~isequal(size(rxns{i}),[1,1])
        multirxns=rxns{i};
        for k=1:length(multirxns)
            if sum(contains(modelTemp.rxns, multirxns{k}))==1
                oldgrstring=modelTemp.grRules{find(contains(modelTemp.rxns, multirxns{k}))};
                grtoadd=A{i};
                grrelation=' or ';
                grRulesstring=[oldgrstring, grrelation, grtoadd];
                modelTemp=changeGeneAssoc(modelTemp, multirxns{k}, grRulesstring, true);
            end
        end
    else
        if sum(contains(modelTemp.rxns, rxns{i}))==1
                oldgrstring=modelTemp.grRules{find(contains(modelTemp.rxns, rxns{i}))};
                grtoadd=A{i};
                grrelation=' or ';
                grRulesstring=[oldgrstring, grrelation, grtoadd];
                modelTemp=changeGeneAssoc(modelTemp, rxns{i}, grRulesstring, true);
        end
    end
end

model1=modelTemp;
