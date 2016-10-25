function [Heter_Network] = InteSTtoNet(DataPath,Heter_Network,DName,CandiMiRNAName)
%==========================================================================
% [Heter_Network]=InteSTtoNet(DName,CandiMiRNAName)
% DName--query disease generic name
% CandiMiRNAName--candidate microRNA list
%==========================================================================
load([DataPath,'miRNA&DiseaseName.mat']);
logicIndex=strcmp(DName,Disease_name);
RowIndex=find(logicIndex,1);
Heter_Network=round((1*10^3)*Heter_Network);
Heter_Network((end+1):(end+2),:)=0;
Heter_Network(end-1,(length(miRNA_name)+RowIndex))=1*10^9;
for i=1:length(CandiMiRNAName)
    logicIndex=strcmp(CandiMiRNAName{i},miRNA_name);
    RowIndex=find(logicIndex,1);
    Heter_Network(end,RowIndex)=1*10^9;
end
addMatrix=zeros(2,2);
Heter_Network(:,(end+1):(end+2))=[Heter_Network(end-1:end,:)';addMatrix];
Heter_Network(end,:)=0;
logicIndex=strcmp(DName,Disease_name);
RowIndex=find(logicIndex,1);
Heter_Network((length(miRNA_name)+RowIndex),end-1)=0;
end
        


