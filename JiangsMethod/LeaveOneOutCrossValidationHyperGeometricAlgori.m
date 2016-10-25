function [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationHyperGeometricAlgori(DataPath,MiFunSim,DiPheSim)
%====================================================================================================
% [CandiRankList,AUC] = LeaveOneOutCrossValidationHyperGeometricAlgori(DataPath,MiFunSim,DiPheSim)
%====================================================================================================

load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
N=length(miRNA_name);
FinalCandiRankList=cell(0,0);
for i=1:length(Disease_name) %#ok<USENS>
    RowArray=find(DiPheSim(i,:));
    miRNAset=cell(0,0);
    for k=1:length(RowArray)
        LogicMiRNArow=strcmp(Disease_name{RowArray(k)},miRNA2disease(:,2)); %#ok<NODEF>
        miRNAset=[miRNAset;miRNA2disease(LogicMiRNArow,1)]; %#ok<AGROW>
    end
    CumProb=zeros(length(miRNA_name),1);
    logicIndex=strcmp(Disease_name{i},miRNA2disease(:,2));
    miRNABin=miRNA2disease(logicIndex,1);
    for j=1:length(miRNABin)
        if j==1
           TempSet=unique([miRNAset;miRNABin((j+1):end)]);
           M=length(TempSet);
        else
           TempSet=unique([miRNAset;miRNABin(1:(j-1));miRNABin((j+1):end)]);
           M=length(TempSet);
        end
        RowIndex=find(strcmp(miRNABin{j},miRNA_name),1);
        MiRNAIndex=find(MiFunSim((RowIndex),:));
        n=length(MiRNAIndex);
        m=length(intersect(TempSet,miRNA_name(MiRNAIndex)));
        Prob=hygepdf(m:n,N,M,n);
        CumProb(j)=1-sum(Prob);
    end
    OtherMiRNA=setdiff(miRNA_name,miRNABin);
    for k=1:length(OtherMiRNA)
        TempSet=unique(miRNAset);
        M=length(TempSet);
        RowIndex=find(strcmp(miRNABin{j},miRNA_name),1);
        MiRNAIndex=find(MiFunSim(RowIndex,:));
        n=length(MiRNAIndex);
        m=length(intersect(TempSet,miRNA_name(MiRNAIndex)));
        Prob=hygepdf(m:n,N,M,n);
        CumProb(length(miRNABin)+k)=1-sum(Prob);
    end
    FinalIndex1=CumProb>0.05;
    FinalIndex2=CumProb(1:length(miRNABin))>0.05;
    FinalIndex3=CumProb((length(miRNABin)+1):end)>0.05;
    if find(FinalIndex1)
        TempCandiRankList=[miRNABin(FinalIndex2);OtherMiRNA(FinalIndex3)];
        TempCandiRankList(:,2)=Disease_name(i);
        TempCandiRankList(:,3)=num2cell(CumProb(FinalIndex1));
        FinalCandiRankList=[FinalCandiRankList;TempCandiRankList];   %#ok<AGROW>
    end
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




