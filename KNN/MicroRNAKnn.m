function [CandiRankList] = MicroRNAKnn(miRNA_name,MiFunSim,mi2diNetwork,DiseaseIndex,K_value)
CandiValueList=zeros(length(miRNA_name),1);
for i=1:length(miRNA_name)
    [~,TopIndex]=sort(MiFunSim(:,i),'descend');
    TopIndex=TopIndex(1:K_value);
    DiseaseMicroRNAIndex=find(mi2diNetwork(:,DiseaseIndex));
    CommonIndex=intersect(TopIndex,DiseaseMicroRNAIndex);
    TempValue=0;
    for j=1:length(CommonIndex)
        TempValue=TempValue+MiFunSim(CommonIndex(j),i)*mi2diNetwork(CommonIndex(j),DiseaseIndex);
    end
    CandiValueList(i)=TempValue;
end
CandiRankList(:,1)=miRNA_name;
CandiRankList(:,3)=num2cell(CandiValueList);
end

