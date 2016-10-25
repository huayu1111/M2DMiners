function [CandiRankListNetCBI] = PredictNetCBI(DataPath,DiseaseName)
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'mi2diNetwork.mat']);
load([DataPath,'DiPheSim.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);
RowIndex=find(strcmp(DiseaseName,Disease_name),1);
AverFunSim=CalAverFunSim(DataPath,mi2diNetwork);
CandiRankListNetCBI(:,1)=miRNA_name;
TempValueList=corrcoef([AverFunSim',DiPheSim(:,RowIndex)]); %#ok<NODEF>
CandiRankListNetCBI(:,3)=num2cell(TempValueList(1:(size(TempValueList,1)-1),size(TempValueList,2)));
NanIndex= isnan(cell2mat(CandiRankListNetCBI(:,3)));
CandiRankListNetCBI(NanIndex,3)={0};
CandiRankListNetCBI(:,2)={DiseaseName};
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2)); %#ok<NODEF>
miRNABin=miRNA2disease(logicIndex,1);
[~,~,FinalIndex1]=intersect(miRNABin,miRNA_name);
CandiRankListNetCBI=CandiRankListNetCBI(FinalIndex1,:);
FinalIndex2=cell2mat(CandiRankListNetCBI(:,3))>0;
CandiRankListNetCBI=CandiRankListNetCBI(FinalIndex2,:);
end

