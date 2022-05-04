% Here we add exchange reactions as well as essential reactions not predicted by the model generator.

% Gap-filling reactions
addedReactions=[];
%CoA synthesis
ADDRXNcoasynt.rxns={'R04231','R03269', 'R02472'};%
ADDRXNcoasynt.mets={{'C00063','C03492','C00097','C00055','C00013','C04352'},{'C04352','C01134','C00011'},{'C00522','C00006','C00966','C00005','C00080'}};%
ADDRXNcoasynt.stoichCoeffs={[-1,-1,-1,1,1,1],[-1,1,1],[-1,-1,1,1,1]};%
ADDRXNcoasynt.lb=[0,0,-1000];%
ADDRXNcoasynt.ub=[1000, 1000, 1000];%
ADDRXNcoasynt.rxnNames={'(R)-4-phosphopantothenate:L-cysteine ligase','N-[(R)-4-Phosphopantothenoyl]-L-cysteine carboxy-lyase','(R)-Pantoate:NADP+ 2-oxidoreductase'};%
model1=addRxns(model1, ADDRXNcoasynt, 1, 'c', true, true);
addedReactions=[addedReactions; ADDRXNcoasynt.rxns'];
model1=changeGeneAssoc(model1,{'R04231','R03269'},'IB49_16075',true);
model1=changeGeneAssoc(model1,'R02472','IB49_13335',true);

clear ADDRXNcoasynt

%Catalase E-value 0, 89% similarity to G.se
ADDRXN.rxns={'R00009'};
ADDRXN.rxnNames={'hydrogen-peroxide:hydrogen-peroxide oxidoreductase'};
ADDRXN.equations={'2 C00027 <=> C00007 + 2 C00001'};
ADDRXN.lb=0;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'R00009','IB49_00430',true);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN


%TCA cycle completion (Sum reaction of isocitrate dehydrogenase). KEGG says
%there should be a proton as product as well but that's not the case when
%MetaCyc metFormulas are used. 
ADDRXNtca.rxns={'R00267'};
ADDRXNtca.mets={'C00311','C00006','C00026','C00011','C00005'};
ADDRXNtca.stoichCoeffs={[-1,-1,1,1,1]};
ADDRXNtca.lb=0;
ADDRXNtca.ub=1000;
ADDRXNtca.rxnNames={'Isocitrate:NADP+ oxidoreductase (decarboxylating)'};
model1=addRxns(model1, ADDRXNtca, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'R00267','IB49_05560',true);
addedReactions=[addedReactions; ADDRXNtca.rxns'];
clear ADDRXNtca

%Shikimate completion
ADDRXNshik.rxns={'R02413'};
ADDRXNshik.equations={'C00493 + C00006 <=> C02637 + C00005 + C00080'};
ADDRXNshik.lb=-1000;
ADDRXNshik.ub=1000;
model1=addRxns(model1, ADDRXNshik, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'R02413','IB49_04485',true);
addedReactions=[addedReactions; ADDRXNshik.rxns'];
clear ADDRXNshik

%FAD biosynthesis
ADDRXNfad.rxns={'R00549', 'R00161'};
ADDRXNfad.equations={'C00002 + C00255 => C00008 + C00061','C00002 + C00061 => C00013 + C00016'};
ADDRXNfad.lb=[0, 0];
ADDRXNfad.ub=[1000, 1000];
model1=addRxns(model1, ADDRXNfad, 1, 'c', false, false);
model1=changeGeneAssoc(model1,ADDRXNfad.rxns,'IB49_16565',true);
addedReactions=[addedReactions; ADDRXNfad.rxns'];
clear ADDRXNfad

% Proline synthesis (spontaneous)
ADDRXN.rxns={'R03314'};
ADDRXN.equations={'C01165 <=> C03912 + C00001'};
ADDRXN.lb=-1000;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

% Leucine synthesis (spontaneous)
ADDRXN.rxns={'R01652'};
ADDRXN.equations={'C00233 + C00011 <=> C04236'};
ADDRXN.lb=-1000;
ADDRXN.ub=0;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

%Thymidylate synthase, EC 2.1.1.45. 
ADDRXN.rxns={'R02101'};
ADDRXN.equations={'C00365 + C00143 <=> C00415 + C00364'};
ADDRXN.lb=0;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
model1=changeGeneAssoc(model1,ADDRXN.rxns,'IB49_00915',true);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

% Lysine biosynthesis
ADDRXNlys.rxns={'R04365'};
ADDRXNlys.equations={'C00091 + C03972 + C00001 <=> C00010 + C04462'};
ADDRXNlys.lb=0;
ADDRXNlys.ub=1000;
model1=addRxns(model1, ADDRXNlys, 1, 'c', true, false);
model1=changeGeneAssoc(model1,'R04365','IB49_15520',true);
addedReactions=[addedReactions; ADDRXNlys.rxns'];
clear ADDRXNlys

ADDRXNhser.rxns={'R01777'};
ADDRXNhser.equations={'C00091 + C00263 <=> C00010 + C01118'};
ADDRXNhser.lb=0;
ADDRXNhser.ub=1000;
model1=addRxns(model1, ADDRXNhser, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'R01777','IB49_00960',true);
addedReactions=[addedReactions; ADDRXNhser.rxns'];
clear ADDRXNhser

% Glycine biosynthesis. From EC 4.1.2.48. Cannot find when blasting, but
% exists according to Cordova et al (2017) 
ADDRXNdtmp.rxns={'R00751'};
ADDRXNdtmp.equations={'C00188 <=> C00037 + C00084'};
ADDRXNdtmp.lb=0;
ADDRXNdtmp.ub=1000;
model1=addRxns(model1, ADDRXNdtmp, 1, 'c', false, false);
addedReactions=[addedReactions; ADDRXNdtmp.rxns'];
clear ADDRXNdtmp

% Acyl-ACP Thioesterase
ADDRXNfa.rxns={'R04014','R08163', 'R08159', 'R01706'};
ADDRXNfa.equations={'C05223 + C00001 <=> C00229 + C02679', 'C04088 + C00001 => C00229 + C01530', 'C05761 + C00001 => C00229 + C06424', 'C05764 + C00001 => C00249 + C00229'};
ADDRXNfa.lb=[0,0,0,0];
ADDRXNfa.ub=[1000,1000,1000,1000];
model1=addRxns(model1, ADDRXNfa, 1, 'c', true, false);
model1=changeGeneAssoc(model1,ADDRXNfa.rxns,'IB49_00680 or IB49_01350 or IB49_18150',true);
addedReactions=[addedReactions; ADDRXNfa.rxns'];
clear ADDRXNfa

% PPP completion
%Adding the spontaneous R02035, and transaldolase (annotated but not picked
%up by homology analysis)
ADDRXNppp.rxns={'R02035', 'R08575'};
ADDRXNppp.equations={'C01236 + C00001 <=> C00345','C05382 + C00118 <=> C00279 + C00085'};
ADDRXNppp.lb=[0,-1000];
ADDRXNppp.ub=[1000, 1000];
model1=addRxns(model1, ADDRXNppp, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'R08575','IB49_09480',true);
addedReactions=[addedReactions; ADDRXNppp.rxns'];
clear ADDRXNppp

% NAD+ Kinase
ADDRXN.rxns={'R00104'};
ADDRXN.rxnNames={'ATP:NAD+ 2´-phosphotransferase'};
ADDRXN.equations={'C00002 + C00003 <=> C00008 + C00006'};
ADDRXN.lb=0;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
model1=changeGeneAssoc(model1,ADDRXN.rxns,'IB49_05820 or IB49_14470',true); 
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

% ATP:thiamine diphosphotransferase. Predicted by KEGG as IB49_16125 
ADDRXN.rxns={'R00619'};
ADDRXN.equations={'C00002 + C00378 <=> C00020 + C00068'};
ADDRXN.lb=0;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
model1=changeGeneAssoc(model1,ADDRXN.rxns,'IB49_16125',true);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

%Pyrophosphate hydrolysis
ADDRXNppi.rxns={'R00004','R00138'};
ADDRXNppi.equations={'C00013 + C00001 <=> 2 C00009 + C00080','C00536 + C00001 <=> C00013 + C00009'};
ADDRXNppi.lb=[0,0];
ADDRXNppi.ub=[1000,1000];
model1=addRxns(model1, ADDRXNppi, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'R00004','IB49_03015 or IB49_07560',true);
model1=changeGeneAssoc(model1,'R00138','IB49_03015 or IB49_07560',true);
addedReactions=[addedReactions; ADDRXNppi.rxns'];
clear ADDRXNppi


% CO2 to Carbonate
ADDRXNbicarb.rxns={'CO2toH2Carbonate', 'H2CarbonatetoCarbonate'};
ADDRXNbicarb.equations={'C00011 + C00001 => C01353', 'C01353 => C00080 + C00288'};
ADDRXNbicarb.lb=[-1000, -1000];
ADDRXNbicarb.ub=[1000, 1000];
model1=addRxns(model1, ADDRXNbicarb, 1, 'c', true, false);
addedReactions=[addedReactions; ADDRXNbicarb.rxns'];
clear ADDRXNbicarb

% Heme synthesis. Change from tRNA-dependent reaction
ADDRXNheme.rxns={'GlutoPreHeme'};
ADDRXNheme.equations={'C00025 + C00002 + C00005 + C00080 => C03741 + C00006 + C00013 + C00020'};
ADDRXNheme.lb=0;
ADDRXNheme.ub=1000;
model1=addRxns(model1, ADDRXNheme, 1, 'c', false, false);
model1=changeGeneAssoc(model1,'GlutoPreHeme','IB49_05115',true);
addedReactions=[addedReactions; ADDRXNheme.rxns'];
clear ADDRXNheme

%Galactose metabolism. UDP-glucose:alpha-D-galactose-1-phosphate uridylyltransferase. 
ADDRXN.rxns={'R00955'};
ADDRXN.equations={'C00029 + C00446 <=> C00103 + C00052'};
ADDRXN.lb=-1000;
ADDRXN.ub=1000;
model1=addRxns(model1, ADDRXN, 1, 'c', false, false);
model1=changeGeneAssoc(model1,ADDRXN.rxns,'IB49_02500',true);
addedReactions=[addedReactions; ADDRXN.rxns'];
clear ADDRXN

% Acetate export
EXRXNac.rxns={'EX_Acetate'};
EXRXNac.equations={'C00033 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

% Lactate export
EXRXNac.rxns={'EX_Lactate'};
EXRXNac.equations={'C00186 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Pyruvate export
EXRXNac.rxns={'EX_Pyruvate'};
EXRXNac.equations={'C00022 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Citrate export
EXRXNac.rxns={'EX_Citrate'};
EXRXNac.equations={'C00158 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Succinate export
EXRXNac.rxns={'EX_Succinate'};
EXRXNac.equations={'C00042 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Alanine export
EXRXNac.rxns={'EX_Alanine'};
EXRXNac.equations={'C00041 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Asparagine export
EXRXNac.rxns={'EX_Asparagine'};
EXRXNac.equations={'C00152 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Aspartate export
EXRXNac.rxns={'EX_Aspartate'};
EXRXNac.equations={'C00049 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Glutamine export
EXRXNac.rxns={'EX_Glutamine'};
EXRXNac.equations={'C00064 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Glutamate export
EXRXNac.rxns={'EX_Glutamate'};
EXRXNac.equations={'C00025 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Glycine export
EXRXNac.rxns={'EX_Glycine'};
EXRXNac.equations={'C00037 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Isoleucine export
EXRXNac.rxns={'EX_Isoleucine'};
EXRXNac.equations={'C00407 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Leucine export
EXRXNac.rxns={'EX_Leucine'};
EXRXNac.equations={'C00123 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Lysine export
EXRXNac.rxns={'EX_Lysine'};
EXRXNac.equations={'C00047 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Methionine export
EXRXNac.rxns={'EX_Methionine'};
EXRXNac.equations={'C00073 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Phenylalanine export
EXRXNac.rxns={'EX_Phenylalanine'};
EXRXNac.equations={'C00079 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Proline export
EXRXNac.rxns={'EX_Proline'};
EXRXNac.equations={'C00148 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Serine export
EXRXNac.rxns={'EX_Serine'};
EXRXNac.equations={'C00065 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Threonine export
EXRXNac.rxns={'EX_Threonine'};
EXRXNac.equations={'C00188 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Tyrosine export
EXRXNac.rxns={'EX_Tyrosine'};
EXRXNac.equations={'C00082 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Valine export
EXRXNac.rxns={'EX_Valine'};
EXRXNac.equations={'C00183 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Arginine export
EXRXNac.rxns={'EX_Arginine'};
EXRXNac.equations={'C00062 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Cysteine export
EXRXNac.rxns={'EX_Cysteine'};
EXRXNac.equations={'C00097 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Histidine export
EXRXNac.rxns={'EX_Histidine'};
EXRXNac.equations={'C00135 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Tryptophan export
EXRXNac.rxns={'EX_Tryptophan'};
EXRXNac.equations={'C00078 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Fumarate export
EXRXNac.rxns={'EX_Fumarate'};
EXRXNac.equations={'C00122 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%Ethanol export
EXRXNac.rxns={'EX_Ethanol'};
EXRXNac.equations={'C00469 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac

%3HB thioesterase
addrxn.rxns={'3HB_thioesterase'};
addrxn.equations={'C03561 + C00001 <=> C01089 + C00010'};
addrxn.lb=0;
addrxn.ub=1000;
model1=addRxns(model1, addrxn, 1, 'c', false, false);
model1=changeGeneAssoc(model1,addrxn.rxns,'IB49_00680 or IB49_01350 or IB49_1815',true);
clear addmets addrxn

%3HB export
EXRXNac.rxns={'EX_3HB'};
EXRXNac.equations={'C01089 <=> '};
EXRXNac.lb=0;
EXRXNac.ub=1000;
model1=addRxns(model1, EXRXNac, 1, 'c', false, false);
clear EXRXNac
%% Adds exchange reactions for the common components in the LC300 medium.
importRxns

