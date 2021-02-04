%this script translates a gene regulatory network with Boolean logic into a
%biochemical reaction network. The script generates a .bngl file which
%encodes these reactions as input to the BioNetGen stochastic simulation
%software. The reason for having a script like this is that it is tedious
%to write the reaction network bngl file manually. The script input is a
%matrix of gene-gene interactions. "0" - no interaction "1" - activating,
%"-1" - repressing, where the column is the source and the row is the
%target. Then, the target gene is considered "active" when at least one
%activator is bound, and no repressors are bound to it.

%this is a made-up matrix of interactions
GeneInteractions= [-1	0	0	1	1	1	0	0	0	0;
	1	1	0	1	1	1	0	0	0	0;
	0	0	1	-1	1	0	0	0	0	0;
	0	0	-1	1	0	0	0	0	0	0;
	0   0	0	0	0	1	0	0	0	0;
	0	1	0	0	0	0	0	0	0	0;
	1	0	0	0	1	1	1	0	0	0;
	1	0	0	0	1	1	0	1	0	0;
	-1	-1	-1	0	-1	-1	-1	-1	0	0;
	1	1	0	0	1	0	1	1	0	0];
    

NumGenes=size(GeneInteractions,1); %number of genes in the network

GeneList={'A','B','C','D','E','F','G','H','I','J'};
RNAList={'RNAFOXC2','RNAGSC','RNAKLF8','RNAmiR_101','RNAmiR_141','RNAmiR_200a','RNAmiR_200b','RNAmiR_200c','RNAmiR_205','RNAmiR_30c'};

%Description of model parameters: h-TF to DNA binding f-TF unbinding, g-
%RNA production rate, k-RNA degradation rate. The parameter values
%are in time units of per-hour
h=10;
f=100;
g2 = 7; %high transcription rate (gene is "on")
g1 = 0.2; %basal transcription rate
g0 = 0; %low transcription rate (gene is "off")
k=1;
B=2; % degree of cooperativity in TF-binding (binding rate is proportional to mRNA^B)

ParNames={'h','f','g0','g1','g2','k','B'};
ParValues=[h,f,g0,g1,g2,k,B];

InitialCopyNumber=5; %this is set somewhat arbitrarily-the initial number of RNAs


%set the Bionetgen file and open for writing
BNGfilename='TestNetwork.bngl';


fid=fopen(BNGfilename,'w+');

%Write the block of parameters and their numeric values
fprintf(fid,'begin model\nbegin parameters \n');
for i = 1:length(ParValues)
    fprintf(fid,'    %s %e \n',ParNames{i},ParValues(i));
end
fprintf(fid,'end parameters \n');
fprintf(fid,'\n');

%Write the block of molecule types, using tilde notation for the different
%possible binding configurations to each promoter 'p1~0~1' means RNA1 can be
%unbound (0) or bound (1)
fprintf(fid,'begin molecule types \n');
for i=1:NumGenes
    GeneName=GeneList{i};
    tempstr=[GeneName '('];
    for j=1:NumGenes-1
        tempstr=[tempstr 'p' num2str(j) '~0~1,'];
    end
    tempstr=[tempstr 'p' num2str(NumGenes) '~0~1)'];
    tempstr=['    ' tempstr ' \n'];
    fprintf(fid,tempstr)
end
%write the RNA types
for i=1:NumGenes
    RNAName=RNAList{i};
    fprintf(fid,['    ' RNAName ' \n']);
end
fprintf(fid,'end molecule types\n');

%Write the block of species, which sets the initial conditions
fprintf(fid,'begin species \n');
%this loop sets all the genes in the fully unbound configuration to start
for i=1:NumGenes
    GeneName=GeneList{i};
    tempstr=[GeneName '('];
    for j=1:NumGenes-1
        tempstr=[tempstr 'p' num2str(j) '~0,'];
    end
    tempstr=[tempstr 'p' num2str(NumGenes) '~0) 1'];
    tempstr=['    ' tempstr ' \n'];
    fprintf(fid,tempstr)
end
%set the initial copy number of the RNA molecules
for i=1:NumGenes
    RNAName=RNAList{i};
    fprintf(fid,['    ' RNAName ' ' num2str(InitialCopyNumber) ' \n']);
end
fprintf(fid,'end species\n');

%this loop sets the observables, which will be written to the trajectory
%file (.gdat) during the simulation. Currently, it only keeps track of the
%RNA species
fprintf(fid,'begin observables \n');
for i=1:NumGenes
    RNAName=RNAList{i};
    fprintf(fid,['    Molecules ' RNAName ' ' RNAName ' \n']);
end
fprintf(fid,'end observables \n');

%this loop writes the functions--rate parameters which are not constant,
%but rather functions of the state of the system
fprintf(fid,'begin functions \n');
for i=1:NumGenes
    RNAName=RNAList{i};
    fprintf(fid,['    hFunc' RNAName '  h*' RNAName '^B\n']);
end
fprintf(fid,'end functions \n');

%begin to write all possible reactions in the network
fprintf(fid,'begin reaction rules \n');

%start with binding/unbinding reactions
for i=1:NumGenes %loop over all the genes
    %this loop figures out which RNAs are regulators of the target
    %gene, then it specifically enumerates the binding reactions for those
    %regulators
    GeneName=GeneList{i}; %the current gene is the target
    RegulatorInds=find(GeneInteractions(i,:)); %find the indices of all genes in the network that regulate the target
    for j=1:numel(RegulatorInds) %loop through the regulators
        RNAName=RNAList{RegulatorInds(j)}; %name the regulator
        tempstr_u='(';
        tempstr_b='(';
        for k=1:NumGenes-1
            if k==RegulatorInds(j)
                tempstr_u=[tempstr_u 'p' num2str(RegulatorInds(j)) '~0,'];
                tempstr_b=[tempstr_b 'p' num2str(RegulatorInds(j)) '~1,'];
            else
                tempstr_u=[tempstr_u 'p' num2str(k) ','];
                tempstr_b=[tempstr_b 'p' num2str(k) ','];
            end
        end
        k=NumGenes;
        if k==RegulatorInds(j)
            tempstr_u=[tempstr_u 'p' num2str(RegulatorInds(j)) '~0)'];
            tempstr_b=[tempstr_b 'p' num2str(RegulatorInds(j)) '~1)'];
        else
            tempstr_u=[tempstr_u 'p' num2str(k) ')'];
            tempstr_b=[tempstr_b 'p' num2str(k) ')'];
        end
        fprintf(fid,['    ' GeneName tempstr_u ' <-> ' GeneName tempstr_b ' hFunc' RNAName ',f\n']);
    end
end

%this block writes the RNA expression (rate g) and degradation (rate k)
%reactions. It implements the rule that the gene is ON only when at
%least one activator and no repressors are bound. (So repressors win).
for i=1:NumGenes
    GeneName=GeneList{i};
    RNAName=RNAList{i};
    ActivatorInds=find(GeneInteractions(i,:)==1); %find the indices of all genes in the network that activate the target
    RepressorInds=find(GeneInteractions(i,:)==-1); %find the indices of all genes in the network that repress the target
    if ActivatorInds %only write an expression reaction if the target has at least one activator
        for j=1:numel(ActivatorInds) %loop through the activators
              tempstr='(';
            for k=1:NumGenes-1
                if k==ActivatorInds(j)
                    tempstr=[tempstr 'p' num2str(ActivatorInds(j)) '~1,'];
                elseif ismember(k,RepressorInds)
                    tempstr=[tempstr 'p' num2str(k) '~0,'];
                else
                    tempstr=[tempstr 'p' num2str(k) ','];
                end
            end
            k=NumGenes;
            if k==ActivatorInds(j)
                tempstr=[tempstr 'p' num2str(ActivatorInds(j)) '~1)'];
            elseif ismember(k,RepressorInds)
                tempstr=[tempstr 'p' num2str(k) '~0)'];
            else
                tempstr=[tempstr 'p' num2str(k) ')'];
            end
            fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g2\n'])
        end
    else
    end
    
    
    %%%Prints Basal Expression for each gene (if all activators and all
    %%%repressors are off then gene X is produced at rate g1)
    tempstr='(';
    for k=1:NumGenes-1
        if ismember(k,ActivatorInds)
             tempstr=[tempstr 'p' num2str(k) '~0,'];
        elseif ismember(k,RepressorInds)
             tempstr=[tempstr 'p' num2str(k) '~0,'];
        else
             tempstr=[tempstr 'p' num2str(k) ','];
        end
    end
    k=NumGenes;
    if ismember(k,ActivatorInds)
        tempstr=[tempstr 'p' num2str(ActivatorInds(j)) '~0)'];
    elseif ismember(k,RepressorInds)
        tempstr=[tempstr 'p' num2str(k) '~0)'];
    else
        tempstr=[tempstr 'p' num2str(k) ')'];
    end
    fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g1\n'])
    
    
    
    %%% Prints "base expression" rate (gene X is produced at rate g0
    %%% regardless of the states of the activators/repressors
    tempstr='(';
    for k=1:NumGenes-1
        tempstr=[tempstr 'p' num2str(k) ','];
    end
    k=NumGenes;
    tempstr=[tempstr 'p' num2str(k) ')'];
    fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g0\n'])
    
            
    
    fprintf(fid,['    ' RNAName ' -> 0 k\n']); %RNA degradation reactions
end


fprintf(fid,'end reaction rules \n');
fprintf(fid,'end model\n');
fprintf(fid,'\n');

fprintf(fid,'generate_network({overwrite=>1})\n');
fprintf(fid,'simulate_ssa({t_start=>0,t_end=>40000,n_steps=>20000})\n');


fclose(fid);

















