function [FinalCandiRankList] = PredictHyperGeometricAlgori(DataPath,MiFunSim,DiPheSim,DiseaseName)
%===============================================================================================
% [FinalCandiRankList] = PredictHyperGeometricAlgori(DataPath,MiFunSim,DiPheSim,DiseaseName)
%===============================================================================================
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
N=length(miRNA_name);
FinalCandiRankList=cell(0,0);
RowIndex=find(strcmp(DiseaseName,Disease_name),1); %#ok<USENS>
RowArray=find(DiPheSim(RowIndex,:));
miRNAset=cell(0,0);
for k=1:length(RowArray)
    LogicMiRNArow=strcmp(Disease_name{RowArray(k)},miRNA2disease(:,2)); %#ok<NODEF>
    miRNAset=[miRNAset;miRNA2disease(LogicMiRNArow,1)]; %#ok<AGROW>
end
CumProb=zeros(length(miRNA_name),1);
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2));
miRNABin=miRNA2disease(logicIndex,1);
miRNAset=[miRNAset;miRNABin];
TempSet=unique(miRNAset);
M=length(TempSet);
for k=1:length(miRNA_name)
    MiRNAIndex=find(MiFunSim(k,:));
    n=length(MiRNAIndex);
    m=length(intersect(TempSet,miRNA_name(MiRNAIndex)));
    Prob=hygepdf(m:n,N,M,n);
    if sum(Prob)<1
        CumProb(k)=1-sum(Prob);
    end
end
[TempCandiRankList,Index]=setdiff(miRNA_name,miRNABin);
TempCandiRankList(:,2)={DiseaseName};
TempCandiRankList(:,3)=num2cell(CumProb(Index));
FinalCandiRankList=[FinalCandiRankList;TempCandiRankList];   
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
end

