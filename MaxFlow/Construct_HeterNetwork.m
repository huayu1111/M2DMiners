function [Heter_Network] = Construct_HeterNetwork(MiFunSim,DiPheSim,MiDiRelat,MiFunField,DiPheField)
%==========================================================================
% Construct_BaseNetwork(MiFunSim,DiPheSim,MiName,DiName,MiDiRelat,MiFamily,MiCluster,Alta,Beta)
% Building base pheno-microRNA assciation network
% MiFunSim--microRNA functional similarity adjacency matrix
% DiPheSim--phenotype semantic similarity adjacency matrix
% MiDiRelat--microRNA and disease relationships
% MiFunField--field value for microRNA functional similarity
% DiPheField--field value for disease phenomenon similarity
%==========================================================================

[MiFunSimNum_Row,~]=size(MiFunSim);
[DiPheSim_Row,~]=size(DiPheSim);
[MiRow,MiCol,~]=find(MiFunSim<MiFunField);
[DiPheRow,DiPheCol,~]=find(DiPheSim<DiPheField);
for i=1:length(MiRow)
    MiFunSim(MiRow(i),MiCol(i))=0;
end
for j=1:length(DiPheRow)
    DiPheSim(DiPheRow(j),DiPheCol(j))=0;
end
MiFunSim=sparse(MiFunSim);
DiPheSim=sparse(DiPheSim);
Heter_Network=sparse((MiFunSimNum_Row+DiPheSim_Row),(MiFunSimNum_Row+DiPheSim_Row));
Heter_Network(1:MiFunSimNum_Row,1:MiFunSimNum_Row)=MiFunSim;
Heter_Network((MiFunSimNum_Row+1):(MiFunSimNum_Row+DiPheSim_Row),(MiFunSimNum_Row+1):(MiFunSimNum_Row+DiPheSim_Row))...
=DiPheSim;

% integrate microRNA-disease association into the heterogenous network
Heter_Network((MiFunSimNum_Row+1):(MiFunSimNum_Row+DiPheSim_Row),1:MiFunSimNum_Row)=MiDiRelat';
Heter_Network(1:MiFunSimNum_Row,(MiFunSimNum_Row+1):(MiFunSimNum_Row+DiPheSim_Row))=MiDiRelat;

end

