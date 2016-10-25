function [FinalCandiRankList,AUC] = LeaveOneOutValidation(DataPath,Heter_Network)
%===========================================================================================
% [FinalCandiRankList,AUC] = LeaveOneOutValidation(DataPath,Heter_Network)
% ==========================================================================================
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
%load([DataPath,'Heter_Network.mat'],'Heter_Network');
FinalCandiRankList=cell(0,0);
for i=1:length(Disease_name) %#ok<USENS>
    logicIndex=strcmp(Disease_name{i},miRNA2disease(:,2)); %#ok<NODEF>
    miRNABin=unique(miRNA2disease(logicIndex,1));
    CandiMiRNAs=setdiff(miRNA_name,miRNABin);
    CandiRankList=cell(0,0);
    for j=1:length(miRNABin)
        TempHeter_Network=Heter_Network;
        TempIndex=strcmp(miRNABin{j},miRNA_name);
        RowIndex=find(TempIndex,1);
        TempHeter_Network(RowIndex,(length(miRNA_name)+i))=0;
        TempHeter_Network((length(miRNA_name)+i),RowIndex)=0;
        Heter_Network_ST=InteSTtoNet(DataPath,TempHeter_Network,Disease_name{i},[miRNABin{j};CandiMiRNAs]);
        [CandiRankValue,~,~,~]=ComputMaxFlowRank(Heter_Network_ST);
        if(isempty(CandiRankValue))
            continue;
        else
            TempCandi=miRNA_name(CandiRankValue(:,1));
            TempCandi(:,2)=Disease_name(i);
            TempCandi(:,3)=num2cell(CandiRankValue(:,2));
            CandiRankList=[CandiRankList;TempCandi]; %#ok<AGROW>
        end
    end
    if(isempty(CandiRankList))
        disp('CandiRankList is empty');
        continue;
    else
        [~,~,IndexB]=intersect(miRNABin,CandiRankList(:,1));
        TempCandiRankList=CandiRankList(IndexB,:);
    end
    for k=1:length(CandiMiRNAs)
        LogicIndex=strcmp(CandiMiRNAs{k},CandiRankList(:,1));
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
    


