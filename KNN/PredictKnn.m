function [CandiRankList] = PredictKnn(DataPath,DiseaseName,mi2diNetwork)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'MiFunSim.mat']);
DiseaseIndex=find(strcmp(DiseaseName,Disease_name),1);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2)); %#ok<NODEF>
miRNABin=miRNA2disease(logicIndex,1);
miRNA_name=setdiff(miRNA_name,miRNABin); %#ok<NODEF>
CandiRankList = MicroRNAKNN(miRNA_name,MiFunSim,mi2diNetwork,DiseaseIndex,20);
CandiRankList(:,2)={DiseaseName};
end