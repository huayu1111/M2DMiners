function [mi2mi] = MiRNAnetwork2Cyto(DataPath)
load([DataPath,'miRNA&DiseaseName.mat'],'miRNA_name');
BooleanMiRNANetwork=ConBooleanMiRNA(DataPath);
mi2mi=cell(0,0);
counter=1;
for i=1:length(BooleanMiRNANetwork)
    ColIndex=find(BooleanMiRNANetwork(i,:));
    for j=1:length(ColIndex)
        mi2mi(counter,1)=miRNA_name(i);
        mi2mi(counter,2)=miRNA_name(ColIndex(j));
        counter=counter+1;
    end
end
end

