function [CandiRankListMBSI] = PredictMBSI(DataPath,DiseaseName)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'mi2diNetwork.mat']);
load([DataPath,'DiPheSim.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
RowIndex=find(strcmp(DiseaseName,Disease_name),1);
AverFunSim=CalAverFunSim(DataPath,mi2diNetwork);
CandiRankListMBSI(:,1)=miRNA_name;
TempValueList=corrcoef([AverFunSim',DiPheSim(:,RowIndex)]); %#ok<NODEF>
CandiRankListMBSI(:,3)=num2cell(TempValueList(1:(size(TempValueList,1)-1),size(TempValueList,2)));
NanIndex= isnan(cell2mat(CandiRankListMBSI(:,3)));
CandiRankListMBSI(NanIndex,3)={0};
CandiRankListMBSI(:,2)={DiseaseName};
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2)); %#ok<NODEF>
miRNABin=miRNA2disease(logicIndex,1);
[~,~,FinalIndex1]=intersect(miRNABin,miRNA_name);
CandiRankListMBSI=CandiRankListMBSI(FinalIndex1,:);
FinalIndex2=cell2mat(CandiRankListMBSI(:,3))>0;
CandiRankListMBSI=CandiRankListMBSI(FinalIndex2,:);
end

