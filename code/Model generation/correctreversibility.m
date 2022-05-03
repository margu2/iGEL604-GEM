% Reaction reversibility corrected based on eQuilibrator  
modelTemp=model;
reversibilitydata=readtable("/Users/emilljungqvist/Documents/Geobacillus/MATLAB/python_json/reversibility_data.xlsx");
reverserxn=[];
forwardrxn=[];
bothrxn=[];
reversibilitydata=table2cell(reversibilitydata);

%% Check for dGm values higher than 30 and lower than -30 Joule/mol and save them to a new array
for i=1:length(reversibilitydata)
    if reversibilitydata{i,2}>=30 
        reverserxn=[reverserxn; reversibilitydata(i,:)];
    elseif reversibilitydata{i,2}<=-30
        forwardrxn=[forwardrxn; reversibilitydata(i,:)];
    else
        bothrxn=[bothrxn; reversibilitydata(i,:)];
    end
end

rev=contains(modelTemp.rxns, reverserxn(:,1));
fwd=contains(modelTemp.rxns, forwardrxn(:,1));
both=contains(modelTemp.rxns, bothrxn(:,1));

for j=1:length(modelTemp.rxns)
    if rev(j)>0
        modelTemp=setParam(modelTemp, 'lb', modelTemp.rxns(j), -1000);
        modelTemp=setParam(modelTemp, 'ub', modelTemp.rxns(j), 0);
        modelTemp.rev(j)=0;
    elseif fwd(j)>0
        modelTemp=setParam(modelTemp, 'lb', modelTemp.rxns(j), 0);
        modelTemp=setParam(modelTemp, 'ub', modelTemp.rxns(j), 1000);
        modelTemp.rev(j)=0;
    elseif both(j)>0
        modelTemp=setParam(modelTemp, 'lb', modelTemp.rxns(j), -1000);
        modelTemp=setParam(modelTemp, 'ub', modelTemp.rxns(j), 1000);
        modelTemp.rev(j)=1;
    end
end

        
model1=modelTemp;
    
