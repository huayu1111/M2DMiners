function [BooleanMiRNANetwork] = ConBooleanMiRNA(DataPath)
load([DataPath,'MiFunSim.mat']);
[MiRow,MiCol,~]=find(MiFunSim>=0.3);
BooleanMiRNANetwork=zeros(length(MiFunSim),length(MiFunSim));
for i=1:length(MiRow)
    BooleanMiRNANetwork(MiRow(i),MiCol(i))=1;
end
end

