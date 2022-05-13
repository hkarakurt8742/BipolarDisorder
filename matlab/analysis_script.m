%%%% REPORTER METABOLITE AND PATHWAY ANALYSIS %%%%

clear 
clc

Network1_Loc = 'Mets_Genes.txt';
Network2_Loc = 'Paths_Mets.txt';
Pvalues_Loc = 'Pvalues.txt';

%load to first Nodes-Edges network
[Network_GM.Mets,Network_GM.Genes]=textread(Network1_Loc,'%s%s');
Network_GM.Genes(1)=[]; Network_GM.Mets(1)=[];% Remove headers
%
%load to second Nodes-Edges network
[Network_PM.Paths,Network_PM.Mets]=textread(Network2_Loc,'%s%s');
Network_PM.Paths(1)=[]; Network_PM.Mets(1)=[];
%load to Omic data
[OmicData.Genes,OmicData.PValues,OmicData.FCs]=textread(Pvalues_Loc,'%s%f%f');
clear Network1_Loc Network2_Loc Pvalues_Loc;

Pvalues.Edges=OmicData.Genes;
Pvalues.PValues=OmicData.PValues;
clear OmicData;
[Network_RG.Rxns,Network_RG.Genes]=textread('Rxns_Genes.txt','%s%s');
Network_RG.Rxns(1)=[]; Network_RG.Genes(1)=[];
k = find(~ismember(Network_RG.Genes, Pvalues.Edges));
Network_RG.Rxns(k)=[]; Network_RG.Genes(k)=[];
k = find(~ismember(Pvalues.Edges, Network_RG.Genes));
Pvalues.Edges(k)=[]; Pvalues.PValues(k)=[];
clear k;

UniRxn = unique(Network_RG.Rxns);
UniGenes = cell(1, 2);
UniGenes{1, 1} = '@@@@@@@@@@';
n=0;
for i = 1:length(UniRxn)
k_1 = find(ismember(Network_RG.Rxns, UniRxn(i)));
k_2 = find(ismember(Pvalues.Edges, Network_RG.Genes(k_1)));
[A, B]= sort(Pvalues.PValues(k_2));
    if length(find(ismember(UniGenes( : , 1), Pvalues.Edges(k_2(B(1)))))) > 0
        continue;
    end
n = n +1;
UniGenes(n, 1) = Pvalues.Edges(k_2(B(1)));
UniGenes{n, 2} = Pvalues.PValues(k_2(B(1)));
clear k_1 k_2 A B;
end
clear i n UniRxn;

Pvalues.Edges = UniGenes( : , 1);
Pvalues.PValues = cell2mat(UniGenes( : , 2));
clear UniGenes;

Network.Nodes=Network_GM.Mets;
Network.Edges=Network_GM.Genes;
clear OmicData Network_GM;

[ Full_Result1 Nodes_Pvalues1 ] = ReporterAnalysis( Network, Pvalues ,100000);

Nodes_Pvalues1 = [Nodes_Pvalues1.Nodes, num2cell(Nodes_Pvalues1.Pvalues)];


E_Pvalues.Edges=Nodes_Pvalues1( : , 1);
E_Pvalues.PValues=cell2mat(Nodes_Pvalues1( : , 2));
Network_N_E.Nodes=Network_PM.Paths;
Network_N_E.Edges=Network_PM.Mets;
clear Network_PM ;

[ Full_Result2 Nodes_Pvalues2 ] = ReporterAnalysis( Network_N_E, E_Pvalues , 100000);
clear Network Pvalues ;
Nodes_Pvalues2 = [Nodes_Pvalues2.Nodes, num2cell(Nodes_Pvalues2.Pvalues)];


%%% GENERATION OF IMAT MODELS %%%


clear
clc

load expression_data

model = readCbModel('iMS570_cobra_with_constraints.xml')
model.c(91) = 0;
model.c(95) = 0;
model.c(97) = 0;
model.c(88) = 1;
model.c(45) = 1;
parsedGPR = GPRparser(model);

hlt_data.gene = genes
hlt_data.value = healthy

bp_data.gene = genes
bp_data.value = bipolar


[gene_id_hlt, gene_expr_hlt] = findUsedGenesLevels(model, hlt_data);

[expressionRxns_hlt parsedGPR] = mapExpressionToReactions(model, hlt_data)

exp_data = [healthy , bipolar];
y = quantile(mean(exp_data') , 9);

hlt_model = iMAT(model, expressionRxns_hlt, y(3) , y(6))

[gene_id_bp, gene_expr_bp] = findUsedGenesLevels(model, bp_data);

[expressionRxns_bp parsedGPR] = mapExpressionToReactions(model, bp_data)

bp_model = iMAT(model, expressionRxns_bp, y(3) , y(6))

model_reactions = model.rxns;

ctrl_rxns = hlt_model.rxns;
bp_rxns = bp_model.rxns;

%%% SAMPOLING %%%


ctrl_rxn_conditions = ismember(model_reactions,ctrl_rxns);
bp_rxn_conditions = ismember(model_reactions,bp_rxns);

rxn_conditions = [ctrl_rxn_conditions bp_rxn_conditions];

hlt_flux = optimizeCbModel(hlt_model);
bp_flux = optimizeCbModel(bp_model);

hlt_model2 = hlt_model;
bp_model2 = bp_model;

hlt_model2.ub(40) = hlt_flux.x(40);
hlt_model2.lb(79) = hlt_flux.x(79)*0.95;

bp_model2.ub(39) = bp_flux.x(39);
bp_model2.lb(78) = bp_flux.x(78)*0.95;

hlt_model2.c(40) = 0;
hlt_model2.c(79) = 0;

bp_model2.c(39) = 0;
bp_model2.c(78) = 0;

% MONTE CARLO SAMPLING -HIT AND RUN
options.nFiles = 1;
options.nStepsPerPoint = 50;
options.nPointsPerFile = 10000;
options.nPointsReturned = 1000;
options.nFilesSkipped = 0;
removeLoopSamplesFlag = false; 

[modelSampling_hlt, samples_hlt] = sampleCbModel(hlt_model2,'healthysampling', 'ACHR',options);

original_reactions = model.rxns;
reduced_rxns = modelSampling_hlt.rxns;
reac_list = ismember(original_reactions , reduced_rxns);
sampling_reactions_hlt = zeros(630,1000);
sampling_reactions_hlt(find(reac_list == 1),:) = samples_hlt;
sampling_reactions_hlt(find(reac_list == 0),:) = 0;

[modelSampling_bp, samples_bp] = sampleCbModel(bp_model2,'bipolarsampling', 'ACHR',options);

original_reactions = model.rxns;
reduced_rxns = modelSampling_bp.rxns;
reac_list = ismember(original_reactions , reduced_rxns);
sampling_reactions_bp = zeros(630,1000);
sampling_reactions_bp(find(reac_list == 1),:) = samples_bp;
sampling_reactions_bp(find(reac_list == 0),:) = 0;

sampling_pvalues = zeros(630,1);

for i = 1:630;
    sampling_pvalues(i,:) = ranksum(sampling_reactions_hlt(i,:) , sampling_reactions_bp(i,:));
end

fdr = mafdr(sampling_pvalues , 'BHFDR' , false);

