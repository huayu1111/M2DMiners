function [CandiRankList] = RLSMDADisease(candiDiseases,Disease_name,MiFunSim,DiPheSim,mi2diNetwork,MiRNAIndex)
w=0.1;
CandiValueList=zeros(length(candiDiseases),1);
[MnRow,~] = size(MiFunSim);
MIMatrix = eye(MnRow);
[DnRow,~] = size(DiPheSim);
DIMatrix = eye(DnRow);
Mres = MiFunSim*(MiFunSim+(MiFunSim*MIMatrix))*mi2diNetwork;
Dres = DiPheSim*(DiPheSim+(DiPheSim*DIMatrix))*(mi2diNetwork');
OptRes = w*Mres+ (1-w)*Dres';
for i=1:length(candiDiseases)
    LogicIndex = strcmp(candiDiseases{i},Disease_name);
    matchIndex = find(LogicIndex,1);
    CandiValueList(i) = OptRes(MiRNAIndex,matchIndex);
end
CandiRankList(:,1)=candiDiseases;
CandiRankList(:,3)=num2cell(CandiValueList);
end
