function [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationMBSI(DataPath)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'MiFunSim.mat']);
load([DataPath,'DiPheSim.mat']);
load([DataPath,'mi2diNetwork.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
DiagElement=diag(DiPheSim); %#ok<NODEF>
DiagElement(:)=1;
DiagElement=diag(DiagElement);
DiPheSim=DiPheSim+DiagElement;
FinalCandiRankList=cell(0,0);
for i=1:length(Disease_name) %#ok<USENS>
    logicIndex=strcmp(Disease_name{i},miRNA2disease(:,2)); %#ok<NODEF>
    miRNABin=unique(miRNA2disease(logicIndex,1));
    TempCandiRankList=cell(0,0);
    for r=1:length(miRNABin)
        Tempmi2diNetwork=mi2diNetwork;
        TempIndex=strcmp(miRNABin{r},miRNA_name);
        RowIndex=find(TempIndex,1);
        Tempmi2diNetwork(RowIndex,i)=0;
        AverFunSim=CalAverFunSim(DataPath,Tempmi2diNetwork);
        TempValueList=corrcoef([AverFunSim',DiPheSim(:,i)]);
        TempValueList=TempValueList(1:(size(TempValueList,1)-1),size(TempValueList,2));
        TempCandi=miRNA_name;
        TempCandi(:,2)=Disease_name(i);
        TempCandi(:,3)=num2cell(TempValueList);
        FinalIndex=cell2mat(TempCandi(:,3))>0;
        TempCandi=TempCandi(FinalIndex,:);
        TempCandiRankList=[TempCandiRankList;TempCandi];      %#ok<AGROW>
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


