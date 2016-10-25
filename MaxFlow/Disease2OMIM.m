function [] = Disease2OMIM(DataPath,infile)
[~,~,Di2OMIM]=xlsread([DataPath,infile]);
load([DataPath,'miRNA&DiseaseName.mat'],'Disease_name');
for i=1:length(Disease_name) %#ok<NODEF>
    logicIndex=strcmp(Disease_name{i},Di2OMIM(:,1));
    Temp=Di2OMIM(logicIndex,3);
    if(~isempty(Temp))
        Disease_name(i,2)=Temp(1); %#ok<AGROW>
    else
        Disease_name(i,2)={''}; %#ok<AGROW>
    end
end
xlswrite('OMIM.xlsx',Disease_name);
end

