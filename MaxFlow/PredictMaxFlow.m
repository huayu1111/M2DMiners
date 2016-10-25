function [CandiRankListMaxFlow] = PredictMaxFlow(DataPath,Heter_Network,DiseaseName)
%==========================================================================
% [CandiRankList] = PredictMaxFlow(DataPath,Heter_Network,DiseaseName)
%==========================================================================
load([DataPath,'miRNA&DiseaseName.mat'],'miRNA_name');
load([DataPath,'miRNA&DiseaseRelationship.mat']);
logicIndex=strcmp(DiseaseName,miRNA2disease(:,2)); %#ok<NODEF>
miRNABin=miRNA2disease(logicIndex,1);
miRNA_nameDiff=setdiff(miRNA_name,miRNABin);
Heter_Network_ST=InteSTtoNet(DataPath,Heter_Network,DiseaseName,miRNA_nameDiff);
[CandiRankValue,~,~,~]=ComputMaxFlowRank(Heter_Network_ST);
CandiRankListMaxFlow=miRNA_name(CandiRankValue(:,1));
CandiRankListMaxFlow(:,3)=num2cell(CandiRankValue(:,2));
CandiRankListMaxFlow(:,2)={DiseaseName};
end

