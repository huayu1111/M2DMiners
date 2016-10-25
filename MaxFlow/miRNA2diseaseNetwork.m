function [mi2diNetwork] = miRNA2diseaseNetwork(DataPath,infile,Alta,Beta,Gamma)
%==========================================================================
% [mi2diNetwork] = miRNA2diseaseNetwork(infile)
% ==input parameters
% infile--microRNA-disease association file
% Alta--microRNA family related network weight parameter
% Beta--microRNA cluster related network weight parameter
% ==output results
% mi2diNetwork--microRNA-disease relationship network
%==========================================================================

load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'MiRNAfamily&cluster.mat']);
[~,mi2di]=xlsread([DataPath,infile]);
miRNA_Num=length(unique(mi2di(:,1)));
Disease_Num=length(unique(mi2di(:,2)));
mi2diNetwork=zeros(miRNA_Num,Disease_Num);
for i=1:length(mi2di)
    microRNA_logicIndex=strcmp(mi2di{i,1},miRNA_name);
    disease_logicIndex=strcmp(mi2di{i,2},Disease_name); %#ok<USENS>
    if(nargin==5)
        mi2diNetwork(microRNA_logicIndex,disease_logicIndex)=Gamma;
    else
        mi2diNetwork(microRNA_logicIndex,disease_logicIndex)=1;
    end       
end

% add miRNA family effect weight to microRNA-disease network
miRNA_family=unique(MiFamily(:,1)); %#ok<NODEF>
for i=1:length(miRNA_family)
    logicIndex1=strcmp(miRNA_family{i,1},MiFamily(:,1));
    miRNA_set1=MiFamily(logicIndex1,2);
    for j=1:length(Disease_name)
        logicIndex2=strcmp(Disease_name{j},mi2di(:,2));
        miRNA_set2=unique(mi2di(logicIndex2,1));
        miRNAinterSet=intersect(miRNA_set1,miRNA_set2);
        if isempty(miRNAinterSet)
            continue;
        else
            for k=1:length(miRNAinterSet)
                microRNA_logicIndex=strcmp(miRNAinterSet{k},miRNA_name);
                microRNA_Row=find(microRNA_logicIndex,1);
                if isempty(microRNA_Row)
                    continue;
                else
                    mi2diNetwork(microRNA_Row,j)=mi2diNetwork(microRNA_Row,j)+Gamma*(length(miRNAinterSet)/(Alta*length(miRNA_set1)));
                end
            end
        end
    end
end

% add miRNA cluster effect weight to microRNA-disease network
miRNA_cluster=unique(MiCluster(:,1)); %#ok<NODEF>
for i=1:length(miRNA_cluster)
    logicIndex1=strcmp(miRNA_cluster{i,1},MiCluster(:,1));
    miRNA_set1=MiCluster(logicIndex1,2);
    for j=1:length(Disease_name)
        logicIndex2=strcmp(Disease_name{j},mi2di(:,2));
        miRNA_set2=unique(mi2di(logicIndex2,1));
        miRNAinterSet=intersect(miRNA_set1,miRNA_set2);
        if isempty(miRNAinterSet)
            continue;
        else
            for k=1:length(miRNAinterSet)
                microRNA_logicIndex=strcmp(miRNAinterSet{k},miRNA_name); 
                microRNA_Row=find(microRNA_logicIndex,1);
                if isempty(microRNA_Row)
                    continue;
                else
                    mi2diNetwork(microRNA_Row,j)=mi2diNetwork(microRNA_Row,j)+Gamma*(length(miRNAinterSet)/(Beta*length(miRNA_set1)));
                end
            end
        end
    end
end
