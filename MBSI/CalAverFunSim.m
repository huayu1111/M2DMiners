function [AverFunSim]=CalAverFunSim(DataPath,mi2diNetwork)
load([DataPath,'miRNA&DiseaseName.mat']);
%load([DataPath,'miRNA&DiseaseRelationship.mat']);
load([DataPath,'MiFunSim.mat']);
AverFunSim=zeros(length(miRNA_name),length(Disease_name));
for i=1:length(miRNA_name)
    for j=1:length(Disease_name)
        LogicIndex=find(mi2diNetwork(:,j));
        Size=length(find(LogicIndex));
        FunSum=0;
        for k=1:length(LogicIndex)
            FunSum=FunSum+MiFunSim(i,LogicIndex(k));
        end
        AverFunSim(i,j)=FunSum/Size;
    end
end

