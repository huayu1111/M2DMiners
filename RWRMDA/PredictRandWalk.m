function [CandiRankListRandomWalk] = PredictRandWalk(DataPath,DiseaseName,Rate)
load([DataPath,'MiFunSim.mat']);
load([DataPath,'miRNA&DiseaseName.mat'],'miRNA_name');
load([DataPath,'miRNA&DiseaseRelationship.mat']);
Diag=diag(sum(MiFunSim));
MiFunSim=((pinv(Diag)).^(0.5))*MiFunSim*((pinv(Diag)).^(0.5));
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2)); %#ok<*NODEF>
miRNAbin=miRNA2disease(logicIndex,1);
ProbArray=zeros(length(miRNA_name),1);
for j=1:length(miRNAbin)
    RowIndex=find(strcmp(miRNAbin{j},miRNA_name),1);
    ProbArray(RowIndex,1)=1/length(miRNAbin);
end
RowIndex=find(ProbArray==0);
TempProb=ProbArray;
ProbArray=(1-Rate).*(MiFunSim*TempProb)+Rate.*(TempProb);
while(norm(ProbArray-TempProb)>1*10^-10)
   TempProb=ProbArray;
   ProbArray=(1-Rate).*(MiFunSim*TempProb)+Rate.*(TempProb);
end
CandiRankListRandomWalk(:,1)=miRNA_name(RowIndex);
CandiRankListRandomWalk(:,3)=num2cell(ProbArray(RowIndex));
CandiRankListRandomWalk(:,2)={DiseaseName};
FinalIndex=cell2mat(CandiRankListRandomWalk(:,3))>=0;
CandiRankListRandomWalk=CandiRankListRandomWalk(FinalIndex,:);
end

