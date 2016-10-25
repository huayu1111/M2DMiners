function [FinalCandiRankList,AUC] =LeaveOneOutValidationDiseases50(DataPath,Heter_Network,Num)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
FinalCandiRankList=cell(0,0);
for i=1:length(miRNA_name) %#ok<USENS>
    logicIndex=strcmp(miRNA_name{i},miRNA2disease(:,1)); %#ok<NODEF>
    DiseaseBin=unique(miRNA2disease(logicIndex,2));
    CandiDiseases=setdiff(Disease_name,DiseaseBin);
    CandiRankList=cell(0,0);
    TempHeter_Network=Heter_Network;
    for r=1:length(DiseaseBin)
        TempIndex=strcmp(DiseaseBin{r},Disease_name);
        RowIndex=find(TempIndex,1);
        TempHeter_Network(i,(length(miRNA_name)+RowIndex))=0;
        TempHeter_Network((length(miRNA_name)+RowIndex),i)=0;
    end
    for j=1:length(DiseaseBin)
        RandSeed=randperm(length(CandiDiseases));
        Heter_Network_ST=InteSTtoNetDisease(DataPath,TempHeter_Network,miRNA_name{i},[DiseaseBin{j};CandiDiseases(RandSeed(1:Num))]);
        [CandiRankValue,~,~,~]=ComputMaxFlowRank(Heter_Network_ST);
        if(isempty(CandiRankValue))
            continue;
        else
            TempCandi=Disease_name((CandiRankValue(:,1)-length(miRNA_name)));
            TempCandi(:,2)=miRNA_name(i);
            TempCandi(:,3)=num2cell(CandiRankValue(:,2));
            CandiRankList=[CandiRankList;TempCandi]; %#ok<AGROW>
        end
    end
    if(isempty(CandiRankList))
        disp('CandiRankList is empty');
        continue;
    else
        [~,~,IndexB]=intersect(DiseaseBin,CandiRankList(:,1));
        TempCandiRankList=CandiRankList(IndexB,:);
    end
    for k=1:length(CandiDiseases)
        LogicIndex=strcmp(CandiDiseases{k},CandiRankList(:,1));
        RowIndex=find(LogicIndex,1);
        if(isempty(RowIndex))
            continue;
        else
            FlowValue=mean(cell2mat(CandiRankList(LogicIndex,3)));
            TempCandiRankList=[TempCandiRankList;CandiRankList(RowIndex,:)]; %#ok<AGROW>
            TempCandiRankList(end,3)={FlowValue};
        end
    end
    FinalCandiRankList=[FinalCandiRankList;TempCandiRankList]; %#ok<AGROW>
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

