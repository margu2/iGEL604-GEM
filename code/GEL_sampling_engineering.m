%% Setup Matlab
clear
clc
options.useFastFVA = false;
changeCobraSolver('gurobi','all');  % Requires working install of Gurobi.
biomassOptPercent_reference = 95;  % Gives 75% biomass formation as lower bound
productOptPercent_reference = 95;  % Gives 90% product formation as lower bound (at 75% growth)
biomassOptPercent = 75;  % Gives 75% biomass formation as lower bound
productOptPercent = 90;  % Gives 90% product formation as lower bound (at 75% growth)
%% Model import
modelFileName = 'LC300/Models/iGEL601_experimental_data_update_2.mat';

model = readCbModel(modelFileName);

%Shut down unnatural exports
model = changeRxnBounds(model,'EX_3HB', 0, 'b');
model = changeRxnBounds(model,'EX_Ethanol', 0, 'b');
changeObjective(model, 'Biomass', 1);


%% Set up for 3HB
model_3HB = changeRxnBounds(model,'3HB_thioesterase', 1000, 'u');
model_3HB = changeRxnBounds(model_3HB,'EX_3HB', 1000, 'u');

% 3HB from glucose
model_3HB_glc = setUpGlcmodel(model_3HB);

% 3HB from xylose
model_3HB_xyl = setUpXylmodel(model_3HB);

% Prepare for sampling
model_3HB_glc_for_sampling = prepareForSampling(model_3HB_glc, 'Biomass', 'EX_3HB', biomassOptPercent, productOptPercent);
model_3HB_xyl_for_sampling = prepareForSampling(model_3HB_xyl, 'Biomass', 'EX_3HB', biomassOptPercent, productOptPercent);

%% General models for glc and xyl
model=setParam(model, 'ub', {'EX_Succinate'}, 1000)
model_glc = setUpGlcmodel(model);
model_xyl = setUpXylmodel(model);

%%
% Prepare for substrate uptake optimiziation (reference sampling)
model_glc_for_sampling = prepareForSampling(model_glc, 'Biomass', 'EX_Glucose', biomassOptPercent_reference, productOptPercent_reference);
model_xyl_for_sampling = prepareForSampling(model_xyl, 'Biomass', 'EX_Xylose', biomassOptPercent_reference, productOptPercent_reference);

%% Sampling parameters
options.nStepsPerPoint = 200;  % Default
options.nPointsReturned = 1000;  
options.nWarmupPoints = 5000;  % Default
outputDirectory = 'Sampled_data_engineering';

%% Sampling

% 75% biomass lb, 90% 3HB lb, glucose
[M_glc_3HB, X_glc_3HB] = sampleCbModel(model_3HB_glc_for_sampling, 'temp_sampling_output/glc_3HB', 'ACHR', options);
samplingData_3HB_glc = extractSamplingData(M_glc_3HB, X_glc_3HB);
writeSampleData(samplingData_3HB_glc, outputDirectory, 'samples_3HB_glc.mat');

data=samplingData_3HB_glc;
rxns=data.rxns;
meanfluxes=data.mean;
stdevfluxes=data.stdev;
t = table(rxns, meanfluxes);
k = table(rxns, stdevfluxes);
writetable(t, 'Escher/SamplingFluxes_3HB_glc_mean.txt');
writetable(k, 'Escher/SamplingFluxes_3HB_glc_stdev.txt');
%
% 75% biomass lb, 90% 3HB lb, xylose
[M_xyl_3HB, X_xyl_3HB] = sampleCbModel(model_3HB_xyl_for_sampling, 'temp_sampling_output/xyl_3HB', 'ACHR', options);
samplingData_3HB_xyl = extractSamplingData(M_xyl_3HB, X_xyl_3HB);
writeSampleData(samplingData_3HB_xyl, outputDirectory, 'samples_3HB_xyl.mat');

data=samplingData_3HB_xyl;
rxns=data.rxns;
meanfluxes=data.mean;
stdevfluxes=data.stdev;
t = table(rxns, meanfluxes);
k = table(rxns, stdevfluxes);
writetable(t, 'Escher/SamplingFluxes_3HB_xyl_mean.txt');
writetable(k, 'Escher/SamplingFluxes_3HB_xyl_stdev.txt');
%


%%
% Reference Sampling optimizing sugar uptake for determining Z score
model_glc_for_sampling=changeRxnBounds(model_glc_for_sampling,'R01976',0,'u')
% 95% biomass lb, 95% EX_Glc lb, glucose
[M_glc_Bio, X_glc_Bio] = sampleCbModel(model_glc_for_sampling, 'temp_sampling_output/glc_Bio', 'ACHR', options);

samplingData_Bio_glc = extractSamplingData(M_glc_Bio, X_glc_Bio);
writeSampleData(samplingData_Bio_glc, outputDirectory, 'samples_Bio_glc.mat');

data=samplingData_Bio_glc;
rxns=data.rxns;
meanfluxes=data.mean;
stdevfluxes=data.stdev;
t = table(rxns, meanfluxes);
k = table(rxns, stdevfluxes);
writetable(t, 'Escher/SamplingFluxes_Bio_glc_mean.txt');
writetable(k, 'Escher/SamplingFluxes_Bio_glc_stdev.txt');



% 95% biomass lb, 95% EX_Xyl lb, xylose
model_xyl_for_sampling=changeRxnBounds(model_xyl_for_sampling,'R01976',0,'u')
[M_xyl_Bio, X_xyl_Bio] = sampleCbModel(model_xyl_for_sampling, 'temp_sampling_output/xyl_Bio', 'ACHR', options);
samplingData_Bio_xyl = extractSamplingData(M_xyl_Bio, X_xyl_Bio);
writeSampleData(samplingData_Bio_xyl, outputDirectory, 'samples_Bio_xyl.mat');

data=samplingData_Bio_xyl;
rxns=data.rxns;
meanfluxes=data.mean;
stdevfluxes=data.stdev;
t = table(rxns, meanfluxes);
k = table(rxns, stdevfluxes);
writetable(t, 'Escher/SamplingFluxes_Bio_xyl_mean.txt');
writetable(k, 'Escher/SamplingFluxes_Bio_xyl_stdev.txt');


%% Function definitions
function samplingModel = prepareForSampling(model, rxn1, rxn2, optPercent1, optPercent2)
%PREPAREFORSAMPLING Local function that sets a model up for sampling
    model = changeObjective(model, rxn1, 1);
    blockedRxns = findBlockedReaction(model);  % Finds reactions that can't carry flux
    model = removeRxns(model, blockedRxns);  % Remove recations that can't carry flux
    
    maxRate1 = optimizeCbModel(model).f;  % Finds max growth rate
    model = changeObjective(model, rxn2, 1);  % Sets objective to maximize 3HB excretion

    %Forces at least optPercent1 % flux through rxn1
    model = changeRxnBounds(model, rxn1, maxRate1*optPercent1/100, 'l');
    
    
    % Remove loops. Change input 3 to 'max' if maximizing biomass or product formtion. Set to 'min' if maximizing sugar uptake (negative flux)     
    [minflux,maxflux] = fluxVariability(model,optPercent2,'min',model.rxns,0,0);  % Gets max and min flux for each reactions without loops
    model = changeRxnBounds(model, model.rxns, maxflux, 'u');  %  Sets upper reaction bounds
    model = changeRxnBounds(model, model.rxns, minflux, 'l');  %  Sets lower reaction bounds
    
    samplingModel = model;
end

function glcModel = setUpGlcmodel(model)
% SETUPGLCMODEL Local function that sets constraints for glc growth
    model = changeRxnBounds(model, 'EX_Glucose', -16.05, 'l');
    model = changeRxnBounds(model, 'EX_Xylose', 0, 'b');
    model = changeRxnBounds(model, 'EX_Oxygen', -34.41, 'l');
    model = changeRxnBounds(model, 'R00086', 10, 'b');  % Maintenance
    
    glcModel = model;
end
    
function xylModel = setUpXylmodel(model)
% SETUPGLCMODEL Local function that sets constraints for xyl growth

    model = changeRxnBounds(model, 'EX_Xylose', -13.61, 'l');
    model = changeRxnBounds(model, 'EX_Glucose', 0, 'b');
    model = changeRxnBounds(model, 'EX_Oxygen', -32.7, 'l');
    model = changeRxnBounds(model,'R00086', 10, 'b');  % Maintenance

    xylModel = model;
end

function sampleDataStruct = extractSamplingData(samplingModel, samplePoints)
%SAMPLEDATASTRUCT Local function that packs sample data into a struct and
%calculates statistics

    struct.model = samplingModel;
    struct.rxns = samplingModel.rxns;
    struct.datapoints = samplePoints;
    struct.stdev = std(samplePoints, 0, 2);
    struct.mean = mean(samplePoints, 2);
    struct.median = median(samplePoints, 2);


    sampleDataStruct = struct; 
end

function writeSampleData(sampleStruct, outDir, fileName)
    outfile = fullfile(outDir, fileName)
    
    structToWrite = sampleStruct;
    save(outfile, 'structToWrite')%sampleStructString)

    tjooohoooo = 1;
end

