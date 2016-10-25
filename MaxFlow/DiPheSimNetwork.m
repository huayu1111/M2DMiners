function [DiPheSim] = DiPheSimNetwork(DataPath,infile1,infile2)
%==========================================================================
% [DiPheSim] = DiPheSimNetwork(DataPath,infile1,infile2)
% ==input parameters
% infile1--disease MeSH term information relationship
% infile2--disease semantic information value
% ==output results
% DiPheSim--disease semantic similarity network
%==========================================================================
[~,~,ID_data]=xlsread([DataPath,infile1]);
[~,~,IC_data]=xlsread([DataPath,infile2]);
load([DataPath,'miRNA&DiseaseName.mat'],'Disease_name');
DiPheSim=zeros(length(Disease_name),length(Disease_name));
for i=1:length(Disease_name)
    for j=(i+1):length(Disease_name)
        RightInfo=0;
        DiseaseNameIndex=strcmp(Disease_name(i),IC_data(:,1));
        RightInfo=RightInfo+IC_data{DiseaseNameIndex,2};
        DiseaseChildrenIndex=strcmp(Disease_name(i),ID_data(:,1));
        DiseaseChildrens=ID_data(DiseaseChildrenIndex,2);
        [RightParentInfo,RightParentArray]=ParentInfo(IC_data,ID_data,DiseaseChildrens);
        RightInfo=RightInfo+RightParentInfo;
        LeftInfo=0;
        DiseaseNameIndex=strcmp(Disease_name(j),IC_data(:,1));
        LeftInfo=LeftInfo+IC_data{DiseaseNameIndex,2};
        DiseaseChildrenIndex=strcmp(Disease_name(j),ID_data(:,1));
        DiseaseChildrens=ID_data(DiseaseChildrenIndex,2);
        [LeftParentInfo,LeftParentArray]=ParentInfo(IC_data,ID_data,DiseaseChildrens);
        LeftInfo=LeftInfo+LeftParentInfo;
        CommonParentArray=intersect(RightParentArray,LeftParentArray);
        CommonInfo=0;
        for k=1:length(CommonParentArray)
            LogicValue=strcmp(CommonParentArray(k),IC_data(:,1));
            CommonInfo=CommonInfo+IC_data{LogicValue,2};
        end
        DiPheSim(i,j)=2*CommonInfo/(RightInfo+LeftInfo);
        DiPheSim(j,i)=2*CommonInfo/(RightInfo+LeftInfo);
    end
    
end


end

