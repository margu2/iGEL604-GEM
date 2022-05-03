%% Creates the Biomass reaction from components 

% Oxidative phosphorylation
addmets.mets={'C00080_e', 'Q', 'QH2'};
addmets.metNames={'H_e', 'Q-cytochrome', 'Q-cytochromeH2'};
addmets.compartments={'e','c','c'};
addmets.metFormulas={'H', 'R', 'RH2'};
addmets.metCharges=[1,0,0];
model1=addMets(model1, addmets);
addrxnoxp.rxns={'Complex_I','Complex_II','Complex_III','ATP_Synthase'};
addrxnoxp.equations={'C00004 + Q + 5 C00080 => C00003 + QH2 + 4 C00080_e','C01352 + Q + C00080 => C00016 + QH2 + C00080_e','0.5 C00007 + QH2 + 2 C00080 => 2 C00080_e + Q + C00001', 'C00008 + C00009 + 4 C00080_e => C00002 + 3 C00080 + C00001'};
addrxnoxp.lb=[0, 0, 0, 0];
addrxnoxp.ub=[1000, 1000, 1000, 1000];
model1=addRxns(model1, addrxnoxp, 1, 'c', true, false);
model1=changeRxns(model1,'R02164', 'C00042 + C00016 + C00080 => C00122 + C01352',1,'c',true);
model1=setParam(model1, 'lb', {'R02164'},0);
model1=setParam(model1, 'ub', {'R02164'},1000);
model1.rev(find(contains(model1.rxns, 'R02164')))=0;
model1=changeGeneAssoc(model1,'ATP_Synthase','IB49_09335 and IB49_09340 and IB49_09345 and IB49_09350 and IB49_09355 and IB49_09360 and IB49_09365 and IB49_09370',true);
model1=changeGeneAssoc(model1,'R02164','IB49_05240',true);
model1=changeGeneAssoc(model1,'Complex_II','IB49_05230 and IB49_05235',true);
model1=changeGeneAssoc(model1,'Complex_III','IB49_15690 and IB49_15695 and IB49_15700 and IB49_15705 and IB49_18085 and IB49_18090',true);
model1=changeGeneAssoc(model1,'Complex_I','IB49_09280 and IB49_09330 and IB49_09325 and IB49_09315 and IB49_09310 and IB49_09305 and IB49_09300 and IB49_09295 and IB49_09290 and IB49_09285 and IB49_09280',true);
clear addmets addrxnoxp



% Cell wall synthesis (45% Peptidoglycan & 55% Teichoic acid). Teichoic
% acid consist of 36 glycerol-phosphates, 1 GlcNAc, 1 ManNAc & each
% glycerol phosphate is branched by either H, d-Ala or Glucose. Teichoic
% acid is trasported out of the cell via an ABC-transporter, hence one ATP
% is hydrolysed per produced acid
% ADDRXNcw.rxns={'Teichoic_acid', 'Peptidoglycan', 'CellWall'};
% ADDRXNcw.equations={'C00043 + C01170 + 36 C00513 + 11.67 C00029 + 11.67 C00993 + C00001 + C00002 => 8.754 Teichoic_acid + C00105 + 12.67 C00015 + 36 C00055 + 11.67 C00133 + C00008 + C00009', 'C05898 => C04574 + Peptidoglycan', '0.55 Teichoic_acid + 0.45 Peptidoglycan => Cellwall'};
% ADDRXNcw.lb=[0,0,0];
% ADDRXNcw.ub=[1000, 1000, 1000];
% model1=addRxns(model1, ADDRXNcw, 1, 'c', true, false);
% addedRxns=[addedRxns; ADDRXNcw.rxns'];
% clear ADDRXNcw


%Teichuronic acid CellWall.
addrxn.rxns={'Teichuronic_acid', 'Peptidoglycan', 'CellWall'};
addrxn.equations={'2 C00203 + 2 C00167 => C00105 + 3 C00015 + 3 C00080 + 0.8366 Teichuronic_acid', 'C05898 => C04574 + Peptidoglycan', '0.55 Teichuronic_acid + 0.45 Peptidoglycan => Cellwall'};
addrxn.lb=[0,0,0];
addrxn.ub=[1000,1000,1000];
model1=addRxns(model1, addrxn, 1, 'c', true, false);
addedReactions=[addedReactions; addrxn.rxns'];
clear addrxn

% Protein synthesis, ATP->AMP for aa activation, 2 GTP -> 2 GDP for peptide
% bond creation. 
ADDRXNprot.rxns={'aa_to_prot'};
ADDRXNprot.equations={'1.075 C00041 + 0.4042 C00152 + 0.4037 C00049 + 0.5674 C00064 + 0.5669 C00025 + 0.8787 C00037 + 0.4789 C00407 + 0.7333 C00123 + 0.5813 C00047 + 0.2087 C00073 + 0.2884 C00079 + 0.3865 C00148 + 0.4000 C00065 + 0.4961 C00188 + 0.2455 C00082 + 0.6649 C00183 + 0.4253 C00062 + 0.1317 C00097 + 0.1727 C00135 + 0.1181 C00078 + 9.2275 C00002 + 18.455 C00044 + 18.455 C00001 => 1.146 Protein + 9.2275 C00020 + 9.2275 C00013 + 18.455 C00035 + 18.455 C00009'};
ADDRXNprot.lb=0;
ADDRXNprot.ub=1000;
model1=addRxns(model1, ADDRXNprot, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNprot.rxns];
clear ADDRXNprot

% Cell membrane Synthesis (tetradecanoic acid, palmitic acid,
% octadecanoic acid and PhosphoGlycerol (2 FAs per PG). 1 water generated per FA bound to PG) 
ADDRXNcm.rxns={'CellMembrane'};
ADDRXNcm.equations={'0.06 C06424 + 2.78 C00249 + 1.02 C01530 + 1.93 C00093 => 1.272 Cellmembrane + 3.86 C00001'};
ADDRXNcm.lb=0;
ADDRXNcm.ub=1000;
model1=addRxns(model1, ADDRXNcm, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNcm.rxns'];
clear ADDRXNcm

% RNA synthesis (nucleotides -> rna + pyrophosphate)
ADDRXNrna.rxns={'RNA'};
ADDRXNrna.equations={'0.726 C00002 + 0.781 C00075 + 0.855 C00063 + 0.756 C00044 => RNA + 3.1 C00013'};
ADDRXNrna.lb=0;
ADDRXNrna.ub=1000;
model1=addRxns(model1, ADDRXNrna, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNrna.rxns'];
clear ADDRXNrna
model1=changeGeneAssoc(model1,'RNA','IB49_10785 and IB49_10790 and IB49_10965 and IB49_16070 and IB49_09505',true);


%DNA synthesis (deoxynucleotides -> dna + pyrophosphate)
ADDRXNdna.rxns={'DNA'};
ADDRXNdna.equations={'0.763 C00131 + 0.824 C00459 + 0.903 C00458 + 0.793 C00286 => DNA + 3.3 C00013'};
ADDRXNdna.lb=0;
ADDRXNdna.ub=1000;
model1=addRxns(model1, ADDRXNdna, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNdna.rxns'];
clear ADDRXNdna
model1=changeGeneAssoc(model1,'DNA','IB49_17680 and IB49_05605 and IB49_16520 and IB49_10325 and IB49_10135 and IB49_04370 and IB49_04210 and IB49_08350 and IB49_10040 and IB49_17005 and IB49_16255 and IB49_05360 and IB49_05535 and IB49_11790 and IB49_17005 and IB49_17495',true);


%Glycogen synthesis (ADP-Glucose -> Glycogen + ADP)
ADDRXNgly.rxns={'Glycogen'};
ADDRXNgly.equations={'C00498 => C00008 + 0.162 Glycogen'};
ADDRXNgly.lb=0;
ADDRXNgly.ub=1000;
model1=addRxns(model1, ADDRXNgly, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'RNA','IB49_06300',true);

addedReactions=[addedReactions; ADDRXNgly.rxns'];
clear ADDRXNgly

% NAD, NADP, adenosylmethinonine, FAD, FMN, PyridoxalPhoshphate, CoA, ThPP, Biotin, Folate, Heme 
ADDRXNcof.rxns={'Cofactor_Pool'};
ADDRXNcof.equations={'0.92 C00003 + 0.054 C00006 + 0.025 C00019 + 0.013 C00016 + 0.022 C00061 + 0.04 C00018 + 0.013 C00010 + 0.024 C00068 + 0.40 C00120 + 0.22 C00504 + 0.16 C00032 => CofactorPool'};
ADDRXNcof.lb=0;
ADDRXNcof.ub=1000;
model1=addRxns(model1, ADDRXNcof, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNcof.rxns'];
clear ADDRXNcof

% Biomass Synthesis: (ATP + Water + Peptidoglycan + Protein + Cell membrane + RNA + DNA) 
ADDRXNbiomass.rxns={'Biomass'};
ADDRXNbiomass.equations={'40 C00002 + 40 C00001 + 0.11 Cellwall + 0.51 Protein + 0.04 Cellmembrane + 0.28 RNA + 0.02 DNA + 0.012 Glycogen + 0.028 CofactorPool => Biomass + 40 C00008 + 40 C00009 + 40 C00080'};
ADDRXNbiomass.lb=0;
ADDRXNbiomass.ub=1000;
model1=addRxns(model1, ADDRXNbiomass, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNbiomass.rxns'];
clear ADDRXNbiomass

%Biomass
EXRXNac.rxns={'EX_Biomass'};
EXRXNac.equations={'Biomass <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac