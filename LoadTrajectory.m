clear
close all

MaxCountLands=20; % Maximum axis length

n=0;
AllCounts=[];

Filename=['C:\Users\clark\OneDrive\Documents\Files\Read_Research\Research_Data\BioNetGen_Sim\Varying_Parameters\Iteration_2\PM_1e2_1_1e2\gdat\0.gdat'];
if isfile(Filename)
    % Import the gdat data into trajectories
    traj=importdata(Filename);
    species=traj.data(:,2:end);
    time=traj.data(:,1);
    titles = split(traj.textdata); % Extract gene names into title array
    titles(1)=[]; % Removes first index which contains '#'
    NumGenes=length(titles)-1; % Subtracts 1 because first title is 'time'
    
    % Plot overall species concentrations vs time
    figure(1)
    plot(time,species)
    
    % Begin computing statistics and subplots of each combination
    NKeepSpecies=NumGenes;
    figure(2)
    for ii=2:NKeepSpecies
        for jj=1:ii-1
            n=n+1; % Index of gene combinations (total = numgenes Choose 2)
            SpeciesInds=[ii,jj];
            X=species(50:end,SpeciesInds(1));
            Y=species(50:end,SpeciesInds(2));
            
            CC(n)=corr(X,Y);
            CC_Spear(n)=corr(X,Y,'Type','Spearman');
            MutualInfoList(n)=mutInfo(X,Y);
            SomeoneZero=find(prod([X,Y],2)==0);
            BothZero=find(sum([X,Y],2)==0);
            TotCellsBothZero=numel(BothZero);
            TotCellsOneZero=numel(SomeoneZero)-TotCellsBothZero;
            notAnotB=numel(BothZero);
            AnotB=numel(find(Y==0))-notAnotB;
            BnotA=numel(find(X==0))-notAnotB;
            AandB=size(X,1)-TotCellsOneZero-notAnotB; %# of cells coexpressing both genes
            
            Coex(n)=AandB/(AnotB+BnotA+AandB);
            OR(n)=AandB*notAnotB/(AnotB*BnotA);
            X(X>MaxCountLands)=MaxCountLands;
            Y(Y>MaxCountLands)=MaxCountLands;
            Xedges=[0:MaxCountLands+1];
            Yedges=Xedges;
            [N,Xedges,Yedges] = histcounts2(X,Y,Xedges,Yedges);
            
            % Creates/selects subplot to plot on, and plots the landscape
            subplot(4,7,n)
            imagesc(N)
            set(gca,'YDir','normal') 
            xlabel(titles(ii+1))
            ylabel(titles(jj+1))
            axis square
            
            LandscapeVec=reshape(N,1,numel(N));
            LandscapeVec=LandscapeVec./sum(LandscapeVec);
            LandscapeArray(n,:)=LandscapeVec;
            Entropy=-sum(LandscapeVec.*log(LandscapeVec),'omitnan');
            EntropyList(n)=Entropy;
        end
        AllCounts=[AllCounts;X];
    end
else
end

%end
% figure(2)
% histogram(AllCounts)
% set(gca,'YScale','log')
% figure(3)
% histogram(CC,30)
% LandscapeArray_NoS=LandscapeArray;
% save LandscapeArray_NoS LandscapeArray_NoS
%
% %FeatureArray=[EntropyList(:),CorrCoeffList(:),MutualInfoList(:),CoexpressionList(:)];
% %save FeatureArray FeatureArray




