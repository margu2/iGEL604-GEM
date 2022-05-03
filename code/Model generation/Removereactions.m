%% Removal of reactions that are wrong or contain Generic compounds

%removeweirdrxns

% PDH complex substep removal 
model1=removeRxns(model1, {'R01699','R02569','R07618', 'R00014'}, 'metFlag', false);

% Acetyl-CoA -> Malonyl-CoA subreactions removal (Acetyl-CoA Carboxylase)
model1=removeRxns(model1, {'R04385', 'R04386'}, 'metFlag', false);

%Biosynthetic reactions using wrong cofactor (for example NADH instead of NADPH)
model1=removeRxns(model1, {'R00705','R10190','R02353','R04444','R01644','R03406','R01836','R08689','R03600','R00361','R00216','R01773','R00746','R07136','R01246','R08282','R08557','R00094','R01565','R00021','R00936','R00093','R00707','R04198','R03191','R00937','R04966','R01257', 'R00848', 'R10507','R00360','R00236','R03050','R04672','R00243'}, 'metFlag', false);
model1=removeRxns(model1, {'R00134','R08017','R00709','R00859','R00430','R04429','R04724','R10221','R10907','R04961','R04958','R04955','R01248','R00842','R04966','R04969','R03291','R07168', 'R02235'}, 'metFlag', false);
model1=removeRxns(model1, {'R09805','R10161','R01058','R00710','R00713','R05705','R03998','R03544','R07417','R10918','R11768','R01429','R00717','R10718','R01520','R02565','R01745','R01388','R02032'}, 'metFlag', false);


% MalCoa + AcCoa -> Dodecanoyl-ACP, doesnt make sense
model1=removeRxns(model1, {'R11671'}, 'metFlag', false);
model1=removeRxns(model1, {'R01625'}, 'metFlag', false);

% PPP with wrong sugars
model1=removeRxns(model1, {'R01827', 'R01830'}, 'metFlag', false);

%2-Oxoglutarate dehydrogenase complex substep removal
model1=removeRxns(model1, {'R01700','R02570','R07618','R00621'}, 'metFlag', false);

%Glycine cleavage system substep removal. Multi-step reaction: R01221
model1=removeRxns(model1, {'R03425','R04125','R03815'}, 'metFlag', false);

%L-glutamate:ferredoxin oxidoreductase (transaminating) substep removal.
%Multi-step reaction: R00021 
model1=removeRxns(model1, {'R10086'}, 'metFlag', false);

%Acyl-ACP Generic Compound reactions
model1=removeRxns(model1, {'R01403', 'R01404', 'R02768', 'R04473', 'R04864', 'R07325', 'R08940', 'R09380', 'R09381', 'R12086', 'R12240'}, 'metFlag', false);

%Beta Oxidation. Cannot be active at the same time as FA biosynthesis,
%since it causes cofactor circulation (NADPH converted to NADH for free)
model1=removeRxns(model1, {'R01975'}, 'metFlag', false);

%Quinone reduction using carbon monoxide. Not likely to happen 
model1=removeRxns(model1, {'R11168'}, 'metFlag', false);

% Reaction using generic compound "Acceptor"
model1=removeRxns(model1, {'R03599'}, 'metFlag', false);

% Reaction removed from KEGG
model1=removeRxns(model1, {'R12174'}, 'metFlag', false);

% Reaction performed in 2 substeps (R04426 & R01652)
model1=removeRxns(model1, {'R10052'}, 'metFlag', false);

% Pyrophosphate-dependent phosphofructokinase. When blasted the sequence
% maps to ATP-dependent pfk of G. Stearothermophilus (E 0.0, Identity 100%), which is the more
% likely candidate. 
model1=removeRxns(model1, {'R00764', 'R02073'}, 'metFlag', false);

% ATP:Carbamate phosphotransferase subreaction removal. Cover reaction
% R00150
model1=removeRxns(model1, {'R01395', 'R12185'}, 'metFlag', false);

% L-Proline:quinone oxidoreductase. Has NADPH-dependant version (R01251)
model1=removeRxns(model1, {'R01253'}, 'metFlag', false);

% Propanoate:CoA ligase (AMP-forming) subreaction removal. COver reaction
% id: R00925
model1=removeRxns(model1, {'R00926', 'R01354'}, 'metFlag', false);

% Succinate:NAD+ oxidoreductase. Succ->Fum conversion carried out by
% electron transport chain reaction R02164
model1=removeRxns(model1, {'R00402'}, 'metFlag', false);

% Homoisocitrate dehydrogenase sub reaction removel. Cover reaction: R01934
model1=removeRxns(model1, {'R04862','R01936'}, 'metFlag', false);

% Dihydrofolate reductase. Cover reaction removal
model1=removeRxns(model1, {'R00940'}, 'metFlag', false);

% (R)-3-Hydroxybutanoate:NAD+ oxidoreductase
model1=removeRxns(model1, {'R01361'}, 'metFlag', false);

%Alpha/Beta-versions of sugars removed.
model1=removeRxns(model1,{'R00802','R00959','R01070','R01600','R01786',...
                          'R01819','R02071','R02703','R02736','R02739',...
                          'R02740','R03321','R04779','R04780','R05133',...
                          'R05134','R05607','R09031','R09780'}, 'metFlag', false);
model1=removeMets(model1, {'C01172','C00668','C05345','C05378','C02336'},false,true);

%Succinyl-CoA -> Succinate transcriptionally regulated. Only active in amino acid degradation 
model1=removeRxns(model1, {'R00410'}, 'metFlag', false)

%Acetate CoA-transferase promiscuous activity causing loops
model1=removeRxns(model1, {'R01179','R01359'}, 'metFlag', false)

%Citrate hydroxymutase. Sub reactions through aconitate used instead. 
model1=removeRxns(model1, {'R01324'}, 'metFlag', false)

%Alanine synthesis using taurine. Other alanine synthesis routes have been
%demonstrated with 13C.
model1=removeRxns(model1, {'R05652'}, 'metFlag', false)


