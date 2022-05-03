% This script works in two steps. First, it finds general metabolites such
% as "Acceptor" and deletes them. It also corrects metabolite formulas so
% that they don't contain . for example.

% The second part of the script tries to find all model metabolites in the
% metacyc database and find charges for them. Surprisingly, many regular
% metabolites, for example D-glucose doesn't seem to exist in the MetaCyc
% database, and will have to be edited manually. Ideally, metabolites
% formulas will also be converted to the MetaCyc version, as the KEGG
% formulas are over-protonated. 

load('metaCycMets2.mat')
Metdatabase=readtable("/Users/emilljungqvist/Documents/Geobacillus/MATLAB/LC300/Data/HomemadeMetDatabase.xlsx");
metdata=[];
NotConverted=[];
MultiHits=[];
metstoberemoved=[];
badformulas=[];
excludemets=[];
excludemetsindex=[];
modelTemp=model1;
badrxns=[];
badformulascorrected=[];
multihits=[];

%Find general metabolites (metabolites with n or R in metformula) and save them to badformulas index 
badformulasindex_n=find(contains(modelTemp.metFormulas, 'n'));
badformulasindex_acc=find(contains(modelTemp.metNames, 'cceptor'));
badformulasindex_R=find(contains(modelTemp.metFormulas, 'R'));
badformulasindex=vertcat(badformulasindex_n, badformulasindex_acc, badformulasindex_R);
badformulas=[modelTemp.mets(badformulasindex), modelTemp.metNames(badformulasindex)];
badformulascorrected=badformulas;

%Change metformula for proton from "p+1" to H
modelTemp.metFormulas(find(contains(modelTemp.metFormulas, 'p+1')))={'H'};

%Remove '.' in metformulas where it exists
modelTemp.metFormulas(find(contains(modelTemp.metFormulas, '.')))=strrep(modelTemp.metFormulas(find(contains(modelTemp.metFormulas, '.'))),'.','');

% List metabolites caught as bad by previous code but shouldn't be removed, and remove them from badformulas array   
excludemets={'[acp]','cyl-carrier protein','Q-cytochrome','Q-cytochromeH2', 'Thioredoxin', 'Thioredoxin disulfide'};
for k=1:length(excludemets)
    excludemetsindex=[find(contains(badformulascorrected(:,2),excludemets{k}))];
    badformulascorrected(excludemetsindex,:)=[];
end

% Further removing generic compounds not caught in previous steps. Remove
% kegg glycan metabolites, as they're already in the model as regular
% compounds
badformulascorrected=vertcat(badformulascorrected,{'C17023', 'Sulfur donor';'G00370','Sucrose (G00370)';'G09795','Trehalose 6-phosphate';...
    'G10552','G10552';'G10553','G10553';'G10555','G10555';'G10556','G10556';'G10610','UDP-N-acetyl-D-glucosamine';'G10614','GDP-D-mannose';...
    'G10619','UDP (G10619)';'G10620','GDP (G10620)';'G11291','Phosphatidyl-myo-inositol monomannoside';'G13102','Phosphatidyl-myo-inositol dimannoside';...
    'G13107','Triacylated phosphatidyl-myo-inositol monomannoside';'G13108','Triacylated phosphatidyl-myo-inositol dimannoside';'C20683','Long-chain acyl-[acyl-carrier protein]'});

%Find the reactions each bad formula metabolite is involved in and save to badrxns
%array.
for j=1:length(badformulascorrected)
    metaboliteindex=find(strcmp(modelTemp.mets, badformulascorrected(j)));
    badrxns=[badrxns; modelTemp.rxns(boolean(modelTemp.S(metaboliteindex,:)))];
end 

%Remove obsolete reactions doesn't work in removeMets for some reason, so the reactions needs to be deleted separately.
modelTemp=removeRxns(modelTemp, badrxns, 'metFlag', false); 
modelTemp=removeMets(modelTemp, badformulascorrected(:,1), false, true);


%% For the remaining metabolites, search for them in the Metacyc database and write their id, name, formula and charge to the metdata table.
for i=1:length(modelTemp.mets)
    if sum(strcmp(metaCycMets.keggid,modelTemp.mets(i)))==1 % Try to match KEGG id of model metabolite to MetaCyc database
        index=find(strcmp(metaCycMets.keggid,modelTemp.mets(i))); % Find the MetaCyc database index of the match and save the data of the metabolite in metdata
        metdata=[metdata; modelTemp.mets(i), modelTemp.metNames(i), metaCycMets.mets(index), metaCycMets.metNames(index), metaCycMets.metFormulas(index), metaCycMets.metCharges(index)];
    elseif sum(strcmp(metaCycMets.metNames, modelTemp.metNames(i)))==1 % Try to match metabolite name of model metabolite to MetaCyc database
        index=find(strcmp(metaCycMets.metNames,modelTemp.metNames(i)));
        metdata=[metdata; modelTemp.mets(i), modelTemp.metNames(i), metaCycMets.mets(index), metaCycMets.metNames(index), metaCycMets.metFormulas(index), metaCycMets.metCharges(index)];
    elseif sum(strcmp(metaCycMets.inchis, modelTemp.inchis(i)))==1 % Try to match inchis of model metabolite to MetaCyc database
        index=find(strcmp(metaCycMets.inchis,modelTemp.inchis(i)));
        metdata=[metdata; modelTemp.mets(i), modelTemp.metNames(i), metaCycMets.mets(index), metaCycMets.metNames(index), metaCycMets.metFormulas(index), metaCycMets.metCharges(index)];
    elseif sum(strcmp(metaCycMets.keggid,modelTemp.mets(i)))>1
        index=find(strcmp(metaCycMets.keggid,modelTemp.mets(i)));
        metdata=[metdata; modelTemp.mets(i), modelTemp.metNames(i), metaCycMets.mets(index(1)), metaCycMets.metNames(index(1)), metaCycMets.metFormulas(index(1)), metaCycMets.metCharges(index(1))];
        for k=1:length(index)
            MultiHits=[MultiHits; modelTemp.mets(i), metaCycMets.metNames(index(k))];
        end
    else 
        metdata=[metdata; modelTemp.mets(i), modelTemp.metNames(i), modelTemp.mets(i), modelTemp.metNames(i), {''}, NaN];
        NotConverted=[NotConverted;modelTemp.mets(i),modelTemp.metNames(i)];
    
    end
end

% Further update the metdata table with the manually precured Metdatabase
for i=1:height(Metdatabase)
    if sum(contains(metdata(:,1), table2cell(Metdatabase(i,1))))==1 % Find corresponding metabolite in metdata
        index=find(contains(metdata(:,1), table2cell(Metdatabase(i,1))));
        metdata(index,:)=table2cell(Metdatabase(i,:)); % Update the row of the metabolite match in metdata with the manually precured data
    elseif sum(contains(metdata(:,1), table2cell(Metdatabase(i,1))))>1
        multihits=[multihits; metdata(find(contains(metdata(:,1), table2cell(Metdatabase(i,1)))))];
        
    end
end

metdatastruct.keggmetid=metdata(:,1);
metdatastruct.keggmetname=metdata(:,2);
metdatastruct.metacycmetid=metdata(:,3);
metdatastruct.metacycmetname=metdata(:,4);
metdatastruct.metformulas=cellstr(metdata(:,5));
metdatastruct.metcharges=cell2mat(metdata(:,6));


modelTemp.metFormulas=metdatastruct.metformulas;
modelTemp.metCharges=metdatastruct.metcharges;
model1=modelTemp;