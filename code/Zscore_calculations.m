
clc
clear

Comparisonfluxes=["R00771","R00756","R01068","R01015","R01061","R01512","R01518",...
                  "R00658","R00200","R00209","R00344","R00351",...
                  "R01325","R00267","R08549","R00405","R02164","R01082",...
                  "R00342","R00479","R00472","R00835","R02035","R01528","R01529","R01056",...
                  "R01641", "R08575","R01067","EX_Acetate","EX_CO2","EX_Oxygen","EX_Biomass","R01858","R01138","R02734"];
              
reference_glc=load('Sampled_data_engineering/samples_Bio_glc.mat');
reference_xyl=load('Sampled_data_engineering/samples_Bio_xyl.mat');
%%
% 3HB production on glucose
glc_3HB=load('2022-03-03_Sampling/Sampled_data_engineering/samples_3HB_glc.mat');
glc_3HB_Zscore=CalculateZscore(glc_3HB, reference_glc, Comparisonfluxes);
writeZscoreData(glc_3HB_Zscore, 'Zscore_data', 'glc_3HB_Zscore');

% 3HB production on xylose
xyl_3HB=load('2022-03-03_Sampling/Sampled_data_engineering/samples_3HB_xyl.mat');
xyl_3HB_Zscore=CalculateZscore(xyl_3HB, reference_xyl, Comparisonfluxes);
writeZscoreData(xyl_3HB_Zscore, 'Zscore_data', 'xyl_3HB_Zscore');


%%
function Z_score_output = CalculateZscore(engineered_model, reference_model, ReactionsToCompare)
Zscore=[];
flux_change_fraction=[];
engineered_mean=[];
reference_mean=[];
engineered_stdev=[];
reference_stdev=[];
%function that calculates the Zscores from datastructs containing sampling
%data.

for i=1:length(ReactionsToCompare)
    engineered_mean=[engineered_mean; engineered_model.structToWrite.mean(find(contains(engineered_model.structToWrite.rxns, ReactionsToCompare(i))))];
    reference_mean=[reference_mean; reference_model.structToWrite.mean(find(contains(reference_model.structToWrite.rxns, ReactionsToCompare(i))))];
    engineered_stdev=[engineered_stdev; engineered_model.structToWrite.stdev(find(contains(engineered_model.structToWrite.rxns, ReactionsToCompare(i))))];
    reference_stdev=[reference_stdev; reference_model.structToWrite.stdev(find(contains(reference_model.structToWrite.rxns, ReactionsToCompare(i))))];

    
end

%Add the equivalent fluxes of PEP-Pyr conversion of Reaction R00200,
    %R001858 & R01138 to make a new mean value and stdev for R00200. Note
    %that the PTS flux will also have to be added to this flux to correctly
    %estimate the conversion rate in the case of growth on glucose
R00200=find(contains(ReactionsToCompare, "R00200"));
R01858=find(contains(ReactionsToCompare, "R01858"));
R01138=find(contains(ReactionsToCompare, "R01138"));
engineered_mean(R00200)=engineered_mean(R00200)+engineered_mean(R01858)+engineered_mean(R01138);
reference_mean(R00200)=reference_mean(R00200)+reference_mean(R01858)+reference_mean(R01138);
engineered_stdev(R00200)=sqrt(engineered_stdev(R00200)^2+engineered_stdev(R01858)^2+engineered_stdev(R01138)^2);
reference_stdev(R00200)=sqrt(reference_stdev(R00200)^2+reference_stdev(R01858)^2+reference_stdev(R01138)^2);

% If R00771 (G6P-F6P) is positive, glucose is used as carbon source. In
% that case, the flux of the PTS system (16.05) is added to R00200 as well.
if reference_mean(1)>=0
    engineered_mean(R00200)=engineered_mean(R00200)-16.05;
    reference_mean(R00200)=reference_mean(R00200)-16.05;
end

%Add the flux of R02734 (lysine biosynthesis substep) to that of R00405 (Succinyl-CoA -> Succinate). 
%During sampling either way is applicable and the flux will be distributed betwen the two.

R00405=find(contains(ReactionsToCompare, "R00405"));
R02734=find(contains(ReactionsToCompare, "R02734"));
%R00405 is defined in the reverse direction (Succinate -> Succinyl-CoA)
engineered_mean(R00405)=-engineered_mean(R00405)+engineered_mean(R02734);
reference_mean(R00405)=-reference_mean(R00405)+reference_mean(R02734);
engineered_stdev(R00405)=sqrt(engineered_stdev(R00405)^2+engineered_stdev(R02734)^2);
reference_stdev(R00405)=sqrt(reference_stdev(R00405)^2+reference_stdev(R02734)^2);


for j=1:length(engineered_mean)
    
    % Calculate Z-score for each reaction. If both fluxes are negative, take
% the absolute of them both to avoid reaction defined in the reverse order
% to fuck shit up. (If a reaction flux is -19 in the reference and -20 in
% the engineered case, the flux has increased, but the Zscore will show a
% negative value indicating a decrease in flux since -20 is a lower number
% than -19)
    if engineered_mean(j)<0 && reference_mean(j)<0
        Zscore=[Zscore; (abs(engineered_mean(j)) - abs(reference_mean(j)))/sqrt(engineered_stdev(j)^2 + reference_stdev(j)^2)];
        
    else
        Zscore=[Zscore; (engineered_mean(j) - reference_mean(j))/sqrt(engineered_stdev(j)^2 + reference_stdev(j)^2)];
    

    end
    
    flux_change_fraction=[flux_change_fraction; engineered_mean(j)/reference_mean(j)];
end


struct.rxns=ReactionsToCompare';
struct.engineered_mean=engineered_mean;
struct.reference_mean=reference_mean;
struct.engineered_stdev=engineered_stdev;
struct.reference_stdev=reference_stdev;
struct.Zscore=Zscore;
struct.flux_change=flux_change_fraction;
Z_score_output=struct;
end


function writeZscoreData(ZscoreStruct, outDir, fileName)
% Writes the Z-score struct to a .mat file, and compiles the reaction means
% and the Z-score in a .txt file
    outfile = fullfile(outDir, fileName);
    ListToWrite = table(ZscoreStruct.rxns, ZscoreStruct.flux_change, ZscoreStruct.Zscore, 'VariableNames',{'Reaction','Flux fractional change', 'Z score'});%
    writetable(ListToWrite, outfile);
    save(outfile, 'ZscoreStruct');

end

