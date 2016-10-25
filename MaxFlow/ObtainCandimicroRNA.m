function [microRNAScore] = ObtainCandimicroRNA(CandiRankValue)
load('miRNA&DiseaseName.mat','miRNA_name');
microRNAScore(:,1)=miRNA_name(CandiRankValue(:,1));
microRNAScore(:,2)=CandiRankValue(:,2);
end

