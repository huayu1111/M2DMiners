function [CandiRankList] = RLSMDA(candiMiRNAs,miRNA_name,MiFunSim,DiPheSim,mi2diNetwork,DiseaseIndex)
w=0.9;
CandiValueList=zeros(length(candiMiRNAs),1);
[MnRow,~] = size(MiFunSim);
MIMatrix = eye(MnRow);
[DnRow,~] = size(DiPheSim);
DIMatrix = eye(DnRow);
Mres = MiFunSim*(MiFunSim+(MiFunSim*MIMatrix))*mi2diNetwork;
Dres = DiPheSim*(DiPheSim+(DiPheSim*DIMatrix))*(mi2diNetwork');
OptRes = w*Mres+ (1-w)*Dres';
for i=1:length(candiMiRNAs)
    LogicIndex = strcmp(candiMiRNAs{i},miRNA_name);
    matchIndex = find(LogicIndex,1);
    CandiValueList(i) = OptRes(matchIndex,DiseaseIndex);
end
CandiRankList(:,1)=candiMiRNAs;
CandiRankList(:,3)=num2cell(CandiValueList);
end

