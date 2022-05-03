% %% Correcting reaction reversibility
% 
% % Pyrophosphate-forming reactions. PPi is quickly degraded to 2 Pi, so
% the reverse of these reactions is unfavourable.
model1=setParam(model1, 'lb', {'R01954','R00578','R00189','R00206',...
                               'R01230','R01231','R00344','R00341','R01049',...
                               'R01280','R02473','R01176','R03383','R01357',...
                               'R00087','R00199'}, 0);
                           
model1=setParam(model1, 'ub', {'R00966','R02142'}, 0);                           
                           
model1.rev(find(contains(model1.rxns, {'R01954','R00578','R00189','R00206',...
                               'R01230','R01231','R00344','R00341','R01049',...
                               'R01280','R02473','R01176','R03383','R01357',...
                               'R00087','R00199','R00966'})))=0;

%Formate:Tetrahydrofolate ligase. dG=-25 at 37 according to MetaCyc. Reverse
%reaction catalyzed by Tetrahydrofolate dehydrogenase.
model1=setParam(model1, 'lb', {'R00943'}, 0);
model1.rev(find(contains(model1.rxns, 'R00943')))=0;


% 
% Amino acid degradation reactions (deaminating). Not irreversible
% in equilibrator but transcriptionally regulated (probably)
model1=setParam(model1, 'lb', {'R01212','R01088', 'R01374', 'R04941', 'R04930', 'R01434','R00490'},0);
model1.rev(find(contains(model1.rxns, 'R01088')))=0;
model1.rev(find(contains(model1.rxns, 'R01374')))=0;
model1.rev(find(contains(model1.rxns, 'R01212')))=0;
model1.rev(find(contains(model1.rxns, 'R04941')))=0;
model1.rev(find(contains(model1.rxns, 'R04930')))=0;
model1.rev(find(contains(model1.rxns, 'R01434')))=0;
model1.rev(find(contains(model1.rxns, 'R00490')))=0;

% Charging of ATP with sulfate for sulfate metabolism. Not reversible
% accoding to equilibrator but biology says otherwise. Probably driven by
% the nextcoming reaction
model1=setParam(model1, 'ub', {'R00529'}, 1000);
model1.rev(find(contains(model1.rxns, 'R00529')))=1;

%D-glucosamine 1,6-phosphomutase. dG=-2.8 at 37C, dG=-29 at 72C
model1=setParam(model1, 'lb', {'R02060'}, -1000);
model1=setParam(model1, 'ub', {'R02060'}, 1000);
model1.rev(find(contains(model1.rxns, 'R02060')))=1;

%UTP:alpha-D-glucose-1-phosphate uridylyltransferase. Pyrophosphate forming
model1=setParam(model1, 'lb', {'R00289'}, 0);
model1=setParam(model1, 'ub', {'R00289'}, 1000);
model1.rev(find(contains(model1.rxns, 'R00289')))=0;

% (S)-dihydroorotate amidohydrolase, dG=3.7 at 37C, dg=-39 at 72C
model1=setParam(model1, 'lb', {'R01993'}, -1000);
model1.rev(find(contains(model1.rxns, 'R01993')))=1;

% Nucleotide synthesis mutase. Irreversible according to equilibrator but
% has to be reversible
model1=setParam(model1, 'ub', {'R07405'}, 1000);
model1.rev(find(contains(model1.rxns, 'R07405')))=1;


% Kinase, has to be driven in the forward direction. Equilibrator is wrong.
% DOI: 10.1074/jbc.274.32.22459
model1=setParam(model1, 'ub', {'R00239','R00480'}, 1000);
model1.rev(find(contains(model1.rxns, 'R00239')))=1;
model1.rev(find(contains(model1.rxns, 'R00480')))=1;

% 
%Glutamine synthase. Regulated by adenylation of the enzyme complex. High
%Glutamine concentrations inactivates the enzyme. doi: 10.1016/j.jmb.2010.09.046
model1=setParam(model1, 'lb', {'R00253'}, 0);
model1.rev(find(contains(model1.rxns, 'R00253')))=0;

%Lysine biosynthesis. dG=-28 at 72C, dG=2.7 at 37C
model1=setParam(model1, 'lb', {'R02291'}, -1000);
model1.rev(find(contains(model1.rxns, 'R02291')))=1;

% Glycine synthase. Chemically reversible but the Glycine cleavage system
% part of the reaction (reverse) is transcriptionally regulated. Only
% active when glycine is at high concentrations
% (https://doi.org/10.1007/BF01659328)
model1=setParam(model1, 'lb', {'R01221'}, 0);
model1.rev(find(contains(model1.rxns, 'R01221')))=0;

% Glutamate dehydrogenase. dG=-39 at 72C, reversible in reality 
model1=setParam(model1, 'lb', {'R00248'}, -1000);
model1.rev(find(contains(model1.rxns, 'R00248')))=1;

% NAD(P) transhydrogenase. NADH->NADPH conversion unfavorable, due to higher concentrations of NADPH 
model1=setParam(model1, 'lb', {'R00112'}, 0);
model1.rev(find(contains(model1.rxns, 'R00112')))=0;

% FADH2:NAD+ oxidoreductase. FADH2->NADH conversion unfavorable, due to higher concentrations of NADH
model1=setParam(model1, 'ub', {'R09748'}, 0);
model1.rev(find(contains(model1.rxns, 'R09748')))=0;

%Malic Enzyme. Reversible in vitro, but irreversible in vivo. Unknown
%regulation
model1=setParam(model1, 'lb', {'R00214'}, 0);
model1.rev(find(contains(model1.rxns, 'R00214')))=0;

%Glyoxylate carboxy-lyase. dG=-20 at 37C, 6 at 72C. Irreversible in aerobic
%conditions
model1=setParam(model1, 'lb', {'R00013'}, 0);
model1.rev(find(contains(model1.rxns, 'R00013')))=0;

%Fructose 1,6 bisphosphatase, only active in gluconeogenesis.
model1=setParam(model1, 'lb', {'R00762','R01845'}, 0);
model1.rev(find(contains(model1.rxns, {'R00762','R01845'})))=0;

%L-glutamate:NADP+ oxidoreductase (transaminating), dG 50 at 37 C. 
model1=setParam(model1, 'ub', {'R00114'}, 0);
model1.rev(find(contains(model1.rxns, 'R00114')))=0;

%L-aspartate ammonia-lyase (fumarate-forming)
%model1=setParam(model1, 'lb', {'R00490'}, 0);
%model1.rev(find(contains(model1.rxns, 'R00490')))=0;

%L-aspartate 1-carboxy-lyase (beta-alanine-forming)
%model1=setParam(model1, 'lb', {'R00489'}, 0);
%model1.rev(find(contains(model1.rxns, 'R00489')))=0;
