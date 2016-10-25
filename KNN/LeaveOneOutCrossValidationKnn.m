function [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationKnn(DataPath,mi2diNetwork)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
load([DataPath,'MiFunSim.mat']);
K_value=20;
FinalCandiRankList=cell(0,0);
for i=1:length(Disease_name) %#ok<USENS>
    logicIndex=strcmp(Disease_name{i},miRNA2disease(:,2)); %#ok<NODEF>
    miRNABin=miRNA2disease(logicIndex,1);
    TempCandiRankList=cell(0,0);
    for j=1:length(miRNABin)
        Tempmi2diNetwork=mi2diNetwork;
        RowIndex=find(strcmp(miRNABin{j},miRNA_name),1);
        Tempmi2diNetwork(RowIndex,i)=0;
        CandiMiRNA=setdiff(miRNA_name,miRNABin);
        TempCandi=MicroRNAKNN([miRNABin{j};CandiMiRNA],MiFunSim,Tempmi2diNetwork,i,K_value);
        TempCandi(:,2)=Disease_name(i);
        TempCandiRankList=[TempCandiRankList;TempCandi]; %#ok<AGROW>
    end
    miRNAUnique=unique(TempCandiRankList(:,1));
    CandiRankList=cell(0,0);
    for s=1:length(miRNAUnique)
        LogicIndex=strcmp(miRNAUnique{s},TempCandiRankList(:,1));
        RowIndex=find(LogicIndex,1);
        RankValue=mean(cell2mat(TempCandiRankList(LogicIndex,3)));
        CandiRankList=[CandiRankList;TempCandiRankList(RowIndex,:)]; %#ok<AGROW>
        CandiRankList(end,3)={RankValue};
    end
    FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>
end
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
GoldData=strcat(miRNA2disease(:,1),miRNA2disease(:,2));
for i=1:length(FinalCandiRankList)
    BridgePair=strcat(FinalCandiRankList{i,1},FinalCandiRankList{i,2});
    if(~isempty(find(strcmp(BridgePair,GoldData),1)));
         FinalCandiRankList(i,4)={1};
         FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    else
         FinalCandiRankList(i,4)={0};
         FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    end
end
AUC=roc_curve(cell2mat(FinalCandiRankList(:,5)),cell2mat(FinalCandiRankList(:,4)));
end

