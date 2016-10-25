function [di2di] = DiPheNetwork2Cyto(DataPath)
load([DataPath,'miRNA&DiseaseName.mat'],'Disease_name');
BooleanPheNetwork=ConBooleanPhe(DataPath);
di2di=cell(0,0);
counter=1;
for i=1:length(BooleanPheNetwork)
    ColIndex=find(BooleanPheNetwork(i,:));
    for j=1:length(ColIndex)
        di2di(counter,1)=miRNA_name(i);
        di2di(counter,2)=miRNA_name(ColIndex(j));
        counter=counter+1;
    end
end
end

