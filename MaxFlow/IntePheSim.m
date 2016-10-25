function [DiPheSim] = IntePheSim(DataPath,infile,DiPheSim)
%================================================================================
% [DiPheSim] = IntePheSim(DataPath,infile,DiPheSim)
%================================================================================
[~,~,Di2OMIM]=xlsread([DataPath,infile]);
PheSimMat=load([DataPath,'MimMiner.mat'],'-ascii');
RowIndex=zeros(1,1);
count=1;
for i=1:size(Di2OMIM,1)
    for j=2:size(Di2OMIM,2)
        if(isnan(cell2mat(Di2OMIM(i,j))))
            break;
        else
            if(~isempty(find(PheSimMat(:,1)==Di2OMIM{i,j},1)))
                RowIndex(count,1)=i;
                RowIndex((count+1),1)=0;
                RowIndex(count,j)=find(PheSimMat(:,1)==Di2OMIM{i,j},1);
            else
                continue;
            end
        end
    end
    if(RowIndex(count,1)~=0)
        count=count+1;
    end
end
RowIndex=RowIndex(1:end-1,:);
cumSum=0;
counter=0;
for k=1:length(RowIndex)
    for s=k+1:length(RowIndex)
        for p=2:length(RowIndex(k,:))
            for q=2:length(RowIndex(s,:))
                if(RowIndex(k,p)~=0&&RowIndex(s,q)~=0)
                    cumSum=cumSum+PheSimMat(RowIndex(k,p),RowIndex(s,q));
                    counter=counter+1;
                end
            end
        end
        SimScore=cumSum/counter;
        if(DiPheSim(RowIndex(k,1),RowIndex(s,1))==0)
            DiPheSim(RowIndex(k,1),RowIndex(s,1))=SimScore;
            DiPheSim(RowIndex(s,1),RowIndex(k,1))=SimScore;
        else
            DiPheSim(RowIndex(k,1),RowIndex(s,1))=(DiPheSim(RowIndex(k,1),RowIndex(s,1))+SimScore)/2;
            DiPheSim(RowIndex(s,1),RowIndex(k,1))=(DiPheSim(RowIndex(k,1),RowIndex(s,1))+SimScore)/2;
        end
        counter=0;
        cumSum=0;
    end
end
end

