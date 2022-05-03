%Removes unused metabolites and genes. 
unusedmets=[];
for i=1:length(model1.mets)
    if sum(abs(model1.S(i,:)))==0
        unusedmets=[unusedmets;model1.mets(i)];
    end
end

%%
model1=removeMets(model1, unusedmets, false, true, true ,true);