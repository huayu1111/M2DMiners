function [MiFunSim] = MiFunSimNetwork(DataPath,infile)
%==========================================================================
% [MiFunSim] = MiFunSimNetwork(infile)
% infile1--microRNA-disease association file
% infile2--microRNA functional similarity file
% MiFunSim--microRNA functional similarity network
%==========================================================================
load([DataPath,'miRNA&DiseaseName.mat'],'miRNA_name');
[~,~,mifunsimi]=xlsread([DataPath,infile]);
MiFunSim=zeros(length(miRNA_name),length(miRNA_name));
for i=1:length(mifunsimi)
        LogicIndex1=strcmp(mifunsimi{i,1},miRNA_name);
        LogicIndex2=strcmp(mifunsimi{i,2},miRNA_name);
        MiFunSim(LogicIndex1,LogicIndex2)=mifunsimi{i,3};
        MiFunSim(LogicIndex2,LogicIndex1)=mifunsimi{i,3};
end
end
