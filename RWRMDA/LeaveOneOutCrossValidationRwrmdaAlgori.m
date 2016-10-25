function [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationRwrmdaAlgori(DataPath,Rate)
load([DataPath,'MiFunSim.mat']);
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
RowSumMiFunSim=sum(MiFunSim);
Diag=diag(RowSumMiFunSim);
MiFunSim=((pinv(Diag)).^(0.5))*MiFunSim*((pinv(Diag)).^(0.5));
FinalCandiRankList=cell(0,0);
for i=1:length(Disease_name) %#ok<*USENS>
    logicIndex=strcmp(Disease_name{i},miRNA2disease(:,2)); %#ok<*NODEF>
    miRNAbin=miRNA2disease(logicIndex,1);
    ProbArray=zeros(length(miRNA_name),1);
    if(length(miRNAbin)==1)
       continue;
    end
    for j=1:length(miRNAbin)
       RowIndex=find(strcmp(miRNAbin{j},miRNA_name),1);
       ProbArray(RowIndex,1)=1/(length(miRNAbin)-1);
    end
    CandiRankList=cell(0,0);
    for k=1:length(miRNAbin)
       TransProbArray=ProbArray;
       Row=find(strcmp(miRNAbin{k},miRNA_name),1);
       TransProbArray(Row,1)=0;
       RowIndex=find(TransProbArray==0);
       TempProb=TransProbArray;
       TransProbArray=(1-Rate).*(MiFunSim*TempProb)+Rate.*(TempProb);
       while(norm(TransProbArray-TempProb)>1*10^-10)
           TempProb=TransProbArray;
           TransProbArray=(1-Rate).*(MiFunSim*TempProb)+Rate.*(TempProb);
       end
       TempCandi=miRNA_name(RowIndex,1);
       TempCandi(:,2)=Disease_name(i);
       TempCandi(:,3)=num2cell(TransProbArray(RowIndex));
       FinalIndex=cell2mat(TempCandi(:,3))>0;
       TempCandi=TempCandi(FinalIndex,:);
       CandiRankList=[CandiRankList;TempCandi]; %#ok<AGROW>
    end
    [~,~,IndexInter]=intersect(miRNAbin,CandiRankList(:,1));
    TempCandiRankList=CandiRankList(IndexInter,:);
    CandiDiff=setdiff(CandiRankList(:,1),miRNAbin);
    for s=1:length(CandiDiff)
        LogicIndex=strcmp(CandiDiff{s},CandiRankList(:,1));
        RowIndex=find(LogicIndex,1);
        RankValue=mean(cell2mat(CandiRankList(LogicIndex,3)));
        TempCandiRankList=[TempCandiRankList;CandiRankList(RowIndex,:)]; %#ok<AGROW>
        TempCandiRankList(end,3)={RankValue};
    end    
    FinalCandiRankList=[FinalCandiRankList;TempCandiRankList]; %#ok<AGROW
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

