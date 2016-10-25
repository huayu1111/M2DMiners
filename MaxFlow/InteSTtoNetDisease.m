function [Heter_Network]=InteSTtoNetDisease(DataPath,Heter_Network,miRNAName,CandiDiseaseName)
load([DataPath,'miRNA&DiseaseName.mat']);
logicIndex=strcmp(miRNAName,miRNA_name);
RowIndex=find(logicIndex,1);
Heter_Network=round((1*10^4)*Heter_Network);
Heter_Network((end+1):(end+2),:)=0;
Heter_Network(end-1,RowIndex)=1*10^9;
for i=1:length(CandiDiseaseName)
    logicIndex=strcmp(CandiDiseaseName{i},Disease_name);
    RowIndex=find(logicIndex,1);
    Heter_Network(end,length(miRNA_name)+RowIndex)=1*10^9;
end
addMatrix=zeros(2,2);
Heter_Network(:,(end+1):(end+2))=[Heter_Network(end-1:end,:)';addMatrix];
Heter_Network(end,:)=0;
logicIndex=strcmp(miRNAName,miRNA_name);
RowIndex=find(logicIndex,1);
Heter_Network(RowIndex,end-1)=0;
end

