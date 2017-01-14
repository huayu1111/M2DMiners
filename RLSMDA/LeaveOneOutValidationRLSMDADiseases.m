function [FinalCandiRankList,AUC] = LeaveOneOutValidationRLSMDADiseases(DataPath,mi2diNetwork)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
load([DataPath,'MiFunSim.mat']);
load([DataPath,'DiPheSim.mat'])
FinalCandiRankList=cell(0,0);
for i=1:length(miRNA_name) %#ok<USENS>
    logicIndex=strcmp(miRNA_name{i},miRNA2disease(:,1)); %#ok<NODEF>
    DiseaseBin=unique(miRNA2disease(logicIndex,2));
    TempCandiRankList=cell(0,0);
    Tempmi2diNetwork = mi2diNetwork;
	for j=1:length(DiseaseBin)
		TempIndex=strcmp(DiseaseBin{j},Disease_name);
		RowIndex=find(TempIndex,1);
		Tempmi2diNetwork(i,RowIndex)=0;
		CandiDiseases=setdiff(Disease_name,DiseaseBin);
		TempCandi=RLSMDADisease([DiseaseBin{j};CandiDiseases],Disease_name,MiFunSim,DiPheSim,Tempmi2diNetwork,i);
		TempCandi(:,2)=miRNA_name(i);
		TempCandiRankList=[TempCandiRankList;TempCandi]; 
	end
    diseaseUnique=unique(TempCandiRankList(:,1));
    CandiRankList=cell(0,0);
    for s=1:length(diseaseUnique)
        LogicIndex=strcmp(diseaseUnique{s},TempCandiRankList(:,1));
        RowIndex=find(LogicIndex,1);
        RankValue=mean(cell2mat(TempCandiRankList(LogicIndex,3)));
        CandiRankList=[CandiRankList;TempCandiRankList(RowIndex,:)]; %#ok<AGROW>
        CandiRankList(end,3)={RankValue};
    end
    FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>
end

% computing the precision, mean rank ratio, concordance and Area under the ROC curve
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
GoldData=strcat(miRNA2disease(:,1),miRNA2disease(:,2));
for i=1:length(FinalCandiRankList)
    BridgePair=strcat(FinalCandiRankList{i,2},FinalCandiRankList{i,1});
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
