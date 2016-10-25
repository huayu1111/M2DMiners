function [mi2diNetwork] =MiRNASetWeight(DataPath,mi2diNetwork,Alta,Beta)
%===========================================================================================
% [mi2diNetwork] =InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,Gamma)
%==========================================================================================
load([DataPath,'miRNA&DiseaseName.mat']);
load([DataPath,'MiRNAfamily&cluster.mat']);
load([DataPath,'miRNA&DiseaseRelationship.mat']);

% add miRNA family effect weight to microRNA-disease network
miRNA_family=unique(MiFamily(:,1)); %#ok<NODEF>
for i=1:length(miRNA_family)
    logicIndex1=strcmp(miRNA_family{i,1},MiFamily(:,1));
    miRNA_set1=MiFamily(logicIndex1,2);
    for j=1:length(Disease_name) %#ok<USENS>
        logicIndex2=strcmp(Disease_name{j},miRNA2disease(:,2)); %#ok<NODEF>
        miRNA_set2=unique(miRNA2disease(logicIndex2,1));
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
                    mi2diNetwork(microRNA_Row,j)=mi2diNetwork(microRNA_Row,j)+length(miRNAinterSet)/(Alta*length(miRNA_set1)); 
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
        logicIndex2=strcmp(Disease_name{j},miRNA2disease(:,2));
        miRNA_set2=unique(miRNA2disease(logicIndex2,1));
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
                    mi2diNetwork(microRNA_Row,j)=mi2diNetwork(microRNA_Row,j)*(1+length(miRNAinterSet)/(Beta*length(miRNA_set1)));
                end
            end
        end
    end
end



