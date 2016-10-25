function [BooleanPheNetwork] = ConBooleanPhe(DataPath)
load([DataPath,'DiPheSim.mat']);
[DiRow,DiCol,~]=find(DiPheSim>=0.3);
BooleanPheNetwork=zeros(length(DiPheSim),length(DiPheSim));
for i=1:length(DiRow)
    BooleanPheNetwork(DiRow(i),DiCol(i))=1;
end
end
