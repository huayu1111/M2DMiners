function [miRNA_family] = FormatTransform(input_file)
%==========================================================================
% [Family] = FormatTransform(input_data) 
% Transform miRNA family or cluster format  
% inputfile--the miRNA family or cluster raw data file
% miRNA_family--the transformed miRNA family
%==========================================================================
dataPath='G:\WorkDir\Data\';
[~,raw_miRNAfamily]=xlsread([dataPath,input_file]);
miRNA_family=cell(0,0);
count=1;
for i=1:length(raw_miRNAfamily)
   row_family=raw_miRNAfamily(i,:);
   for j=2:length(row_family)
           if(~isempty(row_family{1,j}))
                miRNA_family(count,1)=row_family(1,1);
                miRNA_family(count,2)=row_family(1,j);
                count=count+1;
           end
    end     
end

