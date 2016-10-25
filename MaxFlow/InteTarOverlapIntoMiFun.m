function [MiFunSim] = InteTarOverlapIntoMiFun(DataPath,MiFunSim)
%==========================================================================
% [MiFunSim] = InteTarOverlapIntoMiFun(DataPath,MiFunSim)
%==========================================================================
[~,~,MiRNATarOverlap]=xlsread([DataPath,'miRNA_TargetOverlap.xlsx']);
load([DataPath,'miRNA&DiseaseName.mat'],'miRNA_name');
[~,IdxA,IdxB]=intersect(MiRNATarOverlap(:,1),miRNA_name);
MiRNATarOverlap=MiRNATarOverlap(:,2:end);
for i=1:length(IdxB)
    for j=(i+1):length(IdxB)
        if(MiFunSim(IdxB(i),IdxB(j))==0)
            MiFunSim(IdxB(i),IdxB(j))=MiRNATarOverlap{IdxA(i),IdxA(j)};
            MiFunSim(IdxB(j),IdxB(i))=MiRNATarOverlap{IdxA(j),IdxA(i)};
        else
             MiFunSim(IdxB(i),IdxB(j))=(MiRNATarOverlap{IdxA(i),IdxA(j)}+MiFunSim(IdxB(i),IdxB(j)))/2;
             MiFunSim(IdxB(j),IdxB(i))=(MiRNATarOverlap{IdxA(j),IdxA(i)}+MiFunSim(IdxB(j),IdxB(i)))/2;
        end       
    end
end
            
  



