function [ParentInfo,ParentArray] = ParentInfo(IC_data,ID_data,ChildrenNodes)
ParentArray=cell(0,0);
ParentCount=1;
for i=1:length(ChildrenNodes)
    while(length(ChildrenNodes{i})>=3)
        if length(ChildrenNodes{i})==3
            break;
        else
            ChildrenNodes{i}=ChildrenNodes{i}(1:(length(ChildrenNodes{i})-4));
        end
        ChildrenNodeIndex=strcmp(ChildrenNodes{i},ID_data(:,2));
        ParentNodeIndex=strcmp(ID_data(ChildrenNodeIndex,1),IC_data(:,1));
        ParentArray(ParentCount)=IC_data(ParentNodeIndex,1);
        ParentCount=ParentCount+1;      
    end
end
ParentArray=unique(ParentArray);
ParentInfo=0;
for j=1:length(ParentArray)
    InfoIndex=strcmp(ParentArray{j},IC_data(:,1));
    ParentInfo=ParentInfo+IC_data{InfoIndex,2};
end
end

