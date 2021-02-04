function [] = GenerateBNGLModel_fromNetwork(SimN)
%this script translates a gene regulatory network with given Pluripotency genes.
%The script generates a .bngl file which
%encodes these reactions as input to the BioNetGen stochastic simulation
%software. The reason for having a script like this is that it is tedious
%to write the reaction network bngl file manually. The script input is a
%matrix of gene-gene interactions. "0" - no interaction "1" - activating,
%"-1" - repressing, where the column is the source and the row is the
%target. Then, the target gene is considered "active" when at least one
%activator is bound, and no repressors are bound to it. In addition, it
%samples parameters from known parameter distributions 

%SimN = number of files to be generated
NetworkFile=['Nanog_Network' num2str(SimN) '.mat'];
NumGenes=8; %number of genes in the network
Network=zeros(8);
%column: source gene. row: target gene
%order of genes: Gata6,Gcnf,CDX2,KLF4,Nanog,Pbx1,Oct4,Sox2
%Pluripotency network translated from: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006336
Network=[1 0 0 0 -1 0 0 1; 
    1 0 1 0 0 0 0 0;
    0 0 1 0 -1 0 -1 0;
    0 0 0 0 1 0 1 1;
    -1 0 0 1 -1 1 0 0;
    0 0 0 0 1 0 0 0;
    0 -1 -1 0 0 0 1 1;
    0 0 0 0 0 0 1 1];

save(NetworkFile,'Network');

%order of genes: Gata6,Gcnf,CDX2,KLF4,Nanog,Pbx1,Oct4,Sox2
GeneList={'Gata6','Gcnf','Cdx2','Klf4','Nanog','Pbx1','Oct4','Sox2'};
RegList={'gata6','gcnf','cdx2','klf4','nanog','pbx1','oct4','sox2'}; %regulators
RNAList={'rGata6()','rGcnf()','rCdx2()','rKlf4()','rNanog()','rPbx1()','rOct4()','rSox2()'};

%sampling the binding parameters
a=1;
b=5;
f_order=a+(b-a).*rand(1,NumGenes^2);%uniform distribiton of f vals in log-space
f=10.^f_order;

a=log10(f/100);
b=4;
h_order=a+(b-a).*rand(1,NumGenes^2); %randi([-1,2],1,NumGenes);
h=10.^h_order;

save h h
save f f
%sample 'on' expression rates from a log-normal distribution
m=19;
v=350;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

g2 = lognrnd(mu,sigma,1,NumGenes); %high transcription rate (gene is "on")

m=2;v=2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
g1 = lognrnd(mu,sigma,1,NumGenes); %basal transcription rate
g0 = 0; %low transcription rate (gene is "off")
k=1;
C=2; % degree of cooperativity in TF-binding (binding rate is proportional to mRNA^B)


ParNames = [];

%Modifying so that parameters are created for all genes
h_list = [];
f_list = []; 
for i = 1:length(RegList)
    regulator = RegList(i);
    for j = 1:length(GeneList)
        gene = GeneList(j);
        h_list = [h_list, strcat(strcat('h',regulator),gene)];
        f_list = [f_list, strcat(strcat('f',regulator),gene)];
    end
end

ParList = {'g1', 'g2'};
for i = 1:length(ParList)
    for j = 1:NumGenes
        GeneName = GeneList(j);
        ParameterName = strcat(ParList(i),GeneName);
        ParNames = [ParNames, ParameterName];
    end
end
ParNames = [h_list, f_list, ParNames, 'g0', 'k', 'C'];
ParValues=[h,f,g1,g2,g0,k,C];        
        
        
InitialCopyNumber=round(rand(NumGenes,1)*g2(1)); %randomly initializing a subset of genes as 'on'

%set the Bionetgen file and open for writing
BNGfilename=['TestNetwork' num2str(SimN) '.bngl'];

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
% %This writes p01 for every molecule regardless of how many
% activators/repressors

for i=1:NumGenes
    GeneName=GeneList{i};
    tempstr=[GeneName '('];
    for j=1:NumGenes-1
        tempstr=[tempstr 'p' num2str(j) '~0~1,'];
    end
    tempstr=[tempstr 'p' num2str(NumGenes) '~0~1)'];
    tempstr=['    ' tempstr ' \n'];
    fprintf(fid,tempstr);
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
    fprintf(fid,tempstr);
end
%set the initial copy number of the RNA molecules
for i=1:NumGenes
    RNAName=RNAList{i};
    fprintf(fid,['    ' RNAName ' ' num2str(InitialCopyNumber(i)) ' \n']);
end
fprintf(fid,'end species\n');

%this loop sets the observables, which will be written to the trajectory
%file (.gdat) during the simulation. 
fprintf(fid,'begin observables \n');
for i=1:NumGenes
    RNAName=RNAList{i};
    GeneName=GeneList{i};
    fprintf(fid,['    Molecules ' GeneName ' ' RNAName ' \n']);
end
fprintf(fid,'end observables \n');

%this loop writes the functions--rate parameters which are not constant,
%but rather functions of the state of the system
fprintf(fid,'begin functions \n');
% Original code
% for i=1:NumGenes
%     SourceRNAName=RNAList{i};
%     %   SourceGeneName=GeneList{i};
%     SourceRegLett=RegList{i};
%     for j=1:NumGenes
%         %   TargRNAName=RNAList{j};
%         TargGeneName=GeneList{j};
%         fprintf(fid,['    hFunc' SourceRegLett TargGeneName '  h' SourceRegLett TargGeneName '*' SourceRNAName '^C\n']);
%     end
% end

% Edited code to reflect a,b regulator
for i=1:NumGenes
    SourceRNAName=RNAList{i};
    TargGeneName=GeneList{i};
    for j=1:length(RegList)
        SourceRegLett=RegList{j};
        %fprintf(fid,['    hFunc' SourceRegLett TargGeneName '  h*' SourceRegLett TargGeneName '*' SourceRNAName '^C\n']);
        fprintf(fid,['    hFunc' SourceRegLett TargGeneName '  h' SourceRegLett TargGeneName '*' TargGeneName '^C\n']);
    end
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
    RegulatorInds=find(Network(i,:)); %find the indices of all genes in the network that regulate the target
    for j=1:numel(RegulatorInds) %loop through the regulators
        RNAName=RNAList{RegulatorInds(j)}; %name the regulator
        RegLetter=GeneList{RegulatorInds(j)}; %letter of regulator gene
        SourceRegLett=RegList{RegulatorInds(j)};
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
        fprintf(fid,['    ' GeneName tempstr_u ' <-> ' GeneName tempstr_b ' hFunc' SourceRegLett GeneName ',f' SourceRegLett GeneName '\n']);
    end
end

%this block writes the RNA expression (rate g) and degradation (rate k)
%reactions. It implements the rule that the combinatorial expression is:
%if only activators bound, then the transcription rate = N_act*g2
%if only repressors bound, then it is still repressed
%if a combination of repressors and activators, then it is either g2 (more
%activators than repressors), 0 (more repressors than activators), or basal
%(equal number of activators and repressors). If nothing is bound then
%transcription occurs at basal rate

%loop over all genes and their binding states and print the transcription and RNA
%degradation reactions. Note that all possible binding states are printed,
%even though some may never be realized based on allowed binding reactions
%(above)

%CAUTION: The number of nested for-loops depends on the number of genes in network and is
%currently hard-coded!
for i=1:NumGenes %loop over target genes
    GeneName=GeneList{i};
    RNAName=RNAList{i};
%     for Gene1=0:1 %loop over gene 1 bound/unbound to target
%         tempstr1=['(p1' '~' num2str(Gene1) ','];
%         for Gene2=0:1 %loop over gene 2 bound/unbound to target
%             tempstr2=['p2' '~' num2str(Gene2) ')'];
%             tempstr=[tempstr1 tempstr2];
%             BindingState=[Gene1,Gene2];
%             CombinedReg=sum(BindingState.*Network(i,:)); %determine gene regulation state dependent on binding state based on network
%            
%             %applies the correct rna expression equation based on
%             %combinatorial rule
%             if CombinedReg>=1
%                 for ll=1:CombinedReg
%                     fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g2' GeneName '\n']) %activated
%                 end
%             elseif CombinedReg<0
%                 fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g0\n']) %repressed
%             elseif CombinedReg==0
%                 fprintf(fid,['    ' GeneName tempstr ' -> ' GeneName tempstr ' + ' RNAName ' g1' GeneName '\n']) %basal
%             end
%         end
%     end
%     fprintf(fid,['    ' RNAName ' -> 0 k\n']); %RNA degradation reactions
    
    
    currentGeneEffects = Network(i,:); %temp saves the gene effect of current gene in iteration
    indices = find(currentGeneEffects==1|currentGeneEffects==-1); %find indices which hold 1 or -1
    numEffects = length(indices);
    geneStates = zeros(2^numEffects, NumGenes); %predefines array of zeros to hold possible staes: 2^numEffects
    
    %Returns A with all combinations of 1 and 0 (Adapted from
    %online code)
    A = zeros(1,numEffects);
    for k=1:numEffects
        C = nchoosek(1:numEffects,k);
        m = size(C,1);
        R = repmat((1:m)',1,k);
        Ak = accumarray([R(:) C(:)],1,[m numEffects]);
        A = [A; Ak];
    end
    
   %Replaces the predefined geneStates with combinations to form possible
   %states
   for m = 1:(2^numEffects)
       for j = 1:numEffects
           targetIndex = indices(j);
           geneStates(m,targetIndex) = A(m,j);
       end
   end

   %Apply each allowed state to the effect 
   %selects the effect of current gene 
   for k = 1:(2^numEffects)
       %finds the total combined regulating factors
       effect = Network(i,:);
       appliedEffect = effect.*geneStates(k,:);
       CombinedReg = sum(appliedEffect);
       %Writes the promoter string based on 1 or 0 from state
       rule = [];
       for z = 1:NumGenes
           tempstr=['p' num2str(z) '~' num2str(geneStates(k,z)) ','];
           rule = [rule tempstr];
       end
       rule = ['(' rule(1:end-1) ')'];
       %Writes the reaction rules based on the total combined reg fuactors
       if CombinedReg>=1
           for ll=1:CombinedReg
                fprintf(fid,['    ' GeneName rule ' -> ' GeneName rule ' + ' RNAName ' g2' GeneName '\n']); %activated
           end
           elseif CombinedReg<0
                fprintf(fid,['    ' GeneName rule ' -> ' GeneName rule ' + ' RNAName ' g0\n']); %repressed
           elseif CombinedReg==0
                fprintf(fid,['    ' GeneName rule ' -> ' GeneName rule ' + ' RNAName ' g1' GeneName '\n']); %basal
       end
   end
   fprintf(fid,['    ' RNAName ' -> 0 k\n']); %RNA degradation reactions

end

fprintf(fid,'end reaction rules \n');
fprintf(fid,'end model\n');
fprintf(fid,'\n');

fprintf(fid,'generate_network({overwrite=>1})\n');
fprintf(fid,'simulate_ssa({t_start=>0,t_end=>15000,n_steps=>5000})\n');


fclose(fid);
end


    
