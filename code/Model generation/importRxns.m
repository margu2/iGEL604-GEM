importrxns=[];

metstoadd.mets={'C00124_e','C00159_e','C00031_e', 'C00181_e', 'C06232_e','C00853_e','C00314_e','C00378_e','C00255_e','C00253_e','C00568_e','C00120_e','C00504_e','C14818_e','C00009_e','C00093_e','C00095_e'};
metstoadd.metNames={'Galactose_e','Mannose_e','Glucose_e', 'Xylose_e','Molybdate_e','VitB12_e','Pyridoxine_e','Thiamine_e','Ribof_e','Nicotinicacid_e','p-aminobenzoate_e','Biotin_e','Folate_e','Fe_e','Orthophosphate_e','sn-Glycerol 3-phosphate_e','Fructose_e'};
metstoadd.compartments='e';
metstoadd.metFormulas={'C6H12O6','C6H12O6','C6H12O6','C5H10O5','O4Mo','C62H88N13O14PCo','C8H11NO3','C12H17N4OS','C17H20N4O6','C6H4NO2','C7H6NO2','C10H15N2O3S','C19H19N7O6', 'Fe', 'HPO4','C3H7O6P','C6H12O6'};
metstoadd.metCharges=[0,0,0,0,-2,0,0,1,-1,-1,-1,-1,-2,2,-2,-2,0];
model1=addMets(model1, metstoadd);

EXRXNMedium.rxns={'EX_Galactose','EX_Mannose','EX_Xylose', 'EX_Glucose', 'EX_Oxygen', 'EX_Orthophosphate', 'EX_Ammonia', 'EX_Sulfate', 'EX_H2O', 'EX_H', 'EX_CO2'};
EXRXNMedium.equations={'C00124_e <=> ','C00159_e <=> ','C00181_e <=> ', 'C00031_e <=> ', 'C00007 <=> ', 'C00009_e <=> ', 'C00014 <=> ', 'C00059 <=> ', 'C00001 <=> ', 'C00080 <=> ', 'C00011 <=> '};
EXRXNMedium.lb=[-1000,-1000,-1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, 0];
EXRXNMedium.ub=[0,0,0, 0, 0, 1000, 1000, 1000, 1000, 1000, 1000];
model1=addRxns(model1, EXRXNMedium, 1, 'e', false, false);
importrxns=[importrxns; EXRXNMedium.rxns'];
clear EXRXNMedium


EXRXNWolfes.rxns={'EX_Molybdate', 'EX_VitB12', 'EX_PyridoxinHydrochloride', 'EX_Thiamine', 'EX_Nicotinicacid', 'EX_paminobenzoate', 'EX_Biotin', 'EX_Folate', 'EX_Fe'};
EXRXNWolfes.equations={'C06232_e <=> ', 'C00853_e <=> ','C00314_e <=> ','C00378_e <=> ','C00253_e <=> ','C00568_e <=> ','C00120_e <=> ', 'C00504_e <=> ', 'C14818_e <=> '};
EXRXNWolfes.lb=[-1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000];
EXRXNWolfes.ub=[0, 0, 0, 0, 0, 0, 0, 0, 0];
model1=addRxns(model1, EXRXNWolfes, 1, 'e', true, false);
importrxns=[importrxns; EXRXNWolfes.rxns'];
clear EXRXNWolfes

%GLucose import
ADDRXNimport.rxns={'PTS'};
ADDRXNimport.rxnNames={'Phosphotransferase_system'};
ADDRXNimport.equations={'C00022 + C00092 => C00074 + C00031_e'};
ADDRXNimport.lb=-1000;
ADDRXNimport.ub=0;
model1=addRxns(model1, ADDRXNimport, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'PTS','IB49_15245 and IB49_15250 and IB49_09775 and IB49_15255',true);
importrxns=[importrxns; ADDRXNimport.rxns'];
clear ADDRXNimport

%Fructose import
ADDRXNimport.rxns={'PTS_Fru'};
ADDRXNimport.rxnNames={'Phosphotransferase_system_Fructose'};
ADDRXNimport.equations={'C00022 + C01094 => C00074 + C00095_e'};
ADDRXNimport.lb=-30.6;
ADDRXNimport.ub=0;
model1=addRxns(model1, ADDRXNimport, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'PTS_Fru','IB49_01210 and IB49_15255',true);
importrxns=[importrxns; ADDRXNimport.rxns'];
clear ADDRXNimport

%Mannose import. %From blasting of manP of B.subtilis, seems to be
%performed by fructose importer.
ADDRXNimport.rxns={'PTS_Man'};
ADDRXNimport.rxnNames={'Phosphotransferase_system_Mannose'};
ADDRXNimport.equations={'C00022 + C00275 => C00074 + C00159_e'};
ADDRXNimport.lb=-30.6;
ADDRXNimport.ub=0;
model1=addRxns(model1, ADDRXNimport, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'PTS_Man','IB49_01210 and IB49_15255',true);
importrxns=[importrxns; ADDRXNimport.rxns'];
clear ADDRXNimport

%Glycerol import
ADDRXN.rxns={'Glycerol3P_import'};
ADDRXN.equations={'C00093 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00093_e + 2 C00001'};
ADDRXN.lb=-1000;
ADDRXN.ub=0;
model1=addRxns(model1, ADDRXN, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'Glycerol3P_import','IB49_01970 and IB49_01980 and IB49_01985',true);
importrxns=[importrxns; ADDRXN.rxns'];
clear ADDRXN

%ABC importers
ADDRXNimport.rxns={'Xylose_import','Molybdate_import','VitB12_import', 'PyrHCl_import', 'Thiamine_import', 'Nicotinate_import', 'Aminobenzoate_import',...
                   'Biotin_import', 'Folate_import', 'Iron_import', 'Phosphate_import', 'Glucose_import(ABC)','Galactose_import'};
ADDRXNimport.equations={'C00181 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00181_e + 2 C00001','C06232 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C06232_e + 2 C00001',...
                        'C00853 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00853_e + 2 C00001','C00314 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00314_e + 2 C00001',...
                        'C00378 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00378_e + 2 C00001','C00253 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00253_e + 2 C00001',...
                        'C00568 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00568_e + 2 C00001','C00120 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00120_e + 2 C00001',...
                        'C00504 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00504_e + 2 C00001','C14818 + 2 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C14818_e + 2 C00001',...
                        '3 C00009 + 2 C00008 + 2 C00080 => 2 C00002 + C00009_e + 2 C00001', 'C00031 + 2 C00009 + 2 C00008 + 2 C00080 => C00031_e + 2 C00002 + 2 C00001'...
                        'C00124 + 2 C00009 + 2 C00008 + 2 C00080 => C00124_e + 2 C00002 + 2 C00001'};
ADDRXNimport.lb=[-34.9, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -30.6, -1000];
ADDRXNimport.ub=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
model1=addRxns(model1, ADDRXNimport, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'Biotin_import','IB49_06815 or IB49_07190',true);
model1=changeGeneAssoc(model1,'Xylose_import','IB49_01410 and IB49_01390',true);
model1=changeGeneAssoc(model1,'Molybdate_import','IB49_05405 and IB49_14195',true);
model1=changeGeneAssoc(model1,'Phosphate_import','IB49_04110 and IB49_04100 and IB49_04075',true);
model1=changeGeneAssoc(model1,'Glucose_import(ABC)','IB49_08495 and IB49_08490 and IB49_08485 and IB49_13425',true);
%Galactose import has unknown gene relationship. Possibly done by Glucose
%ABC importer or PTS system.
importrxns=[importrxns; ADDRXNimport.rxns'];
clear ADDRXNimport