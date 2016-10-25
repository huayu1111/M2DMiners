function [DiPheSim] = DiPheSimNetwork_old(DataPath,infile1,infile2)
%=============================================================================================== 
% [DiPheSim] = DiPheSimNetwork(DataPath,infile1,infile2)
% building disease phenotype similar network
% infile1--disease MeSH ID information
% infile2--disease information content value
%================================================================================================
[~,~,DiID]=xlsread([DataPath,infile1]);
% DiID(:,(size(DiID,2)+1))={NaN};
% building disease MeSH term tree
DiTree=tree('RootNode');
for i=2:size(DiID,2)
    if(i==2)
        recurID=unique(DiID(:,i));
        logicCell=cellfun(@isnan,recurID,'UniformOutput',false);
        logicLength=cellfun(@length,logicCell);
        logicIndexTotal=find(logicLength~=1);
        relatPair=zeros(length(logicIndexTotal),2);
        for k=1:length(logicIndexTotal)
            logicIndexBin=strcmp(recurID(logicIndexTotal(k)),DiID(:,i));
            DiseaseIndex=find(logicIndexBin,1);
            logicIndexBin=false(length(logicIndexBin),1);
            logicIndexBin(DiseaseIndex)=true;
            % save disease index
            relatPair(k,1)=DiseaseIndex;
            % add disease name to DiTree
            [DiTree,relatPair(k,2)]=DiTree.addnode(1,DiID(logicIndexBin,1));
        end
        ParentNodeID=DiTree.findleaves();
    else
        tempPair=relatPair;
        relatPair=zeros(0,0);
        % save leave nodes
        count=1;
        for j=1:length(ParentNodeID)
            RowIndex=find(tempPair(:,2)==ParentNodeID(j));
            if (i==3)
                LogicIndex=strcmp(DiID{tempPair(RowIndex,1),i-1},DiID(:,i-1));
            else
               TempNodeBin=cellfun(@num2str,DiID(:,i-1),'UniformOutput',false);
               TempNode=num2str(DiID{tempPair(RowIndex,1),i-1});
               LogicIndex=strcmp(TempNode,TempNodeBin);
            end
            LogicIndex(tempPair(RowIndex,1))=false;
            % computing the number of childrens;
            ChildrenNodeID=unique(cellfun(@num2str,DiID(LogicIndex,i),'UniformOutput',false));
            for k=1:length(ChildrenNodeID)
                TempChildrenBin=cellfun(@num2str,DiID(:,i),'UniformOutput',false);
                logicIndexBin=strcmp(ChildrenNodeID(k),TempChildrenBin);
                DiseaseIndex=find(logicIndexBin,1);
                logicIndexBin=false(length(logicIndexBin),1);
                logicIndexBin(DiseaseIndex)=true;
                relatPair(count,1)=DiseaseIndex;
                [DiTree,relatPair(count,2)]=DiTree.addnode(ParentNodeID(j),DiID(logicIndexBin,1));
                count=count+1;
            end
        end
        clear ParentNodeID;
        ParentNodeID=relatPair(:,2);
    end
    fprintf('please waiting, the %d-th layer tree building successful\n',i);
end

% computing disease phenotye similarity based on disease MeSH term information content
[~,~,InfoValue]=xlsread([DataPath,infile2]);
load('miRNA&DiseaseName.mat','Disease_name');
DiPheSim=zeros(length(Disease_name),length(Disease_name));
for i=1:length(Disease_name)
    for j=(i+1):length(Disease_name)
        RightParentID=find(strcmp(DiTree,Disease_name(i)));
        LeftParentID=find(strcmp(DiTree,Disease_name(j)));
        RightParentContent=cell(0,0);
        LeftParentContent=cell(0,0);
        count=1;
        RightInfo=0;
        while RightParentID~=1
            RightParentContent{count,1}=DiTree.get(RightParentID);
            LogicValue=strcmp(RightParentContent{count,1},InfoValue(:,1));
            RightInfo=RightInfo+InfoValue{LogicValue,2};
            RightParentID=find(strcmp(DiTree,RightParentContent{count,1}));
            count=count+1;
        end
        count=1;
        LeftInfo=0;
        while LeftParentID~=1
            LeftParentContent{count,1}=DiTree.get(LeftParentID);
            LogicValue=strcmp(LeftParentContent{count,1},InfoValue(:,1));
            LeftInfo=LeftInfo+InfoValue{LogicValue,2};
            LeftParentID=find(strcmp(DiTree,LeftParentContent{count,1}));
            count=count+1;
        end
        CommonParentContent=intersect(RightParentContent,LeftParentContent);
        CommonInfo=0;
        for k=1:length(CommonParentContent)
            LogicValue=strcmp(CommonParentContent(k),InfoValue(:,1));
            CommonInfo=CommonInfo+InfoValue{LogicValue,2};
        end
        DiPheSim(i,j)=2*CommonInfo/(RightInfo+LeftInfo);
        DiPheSim(j,i)=2*CommonInfo/(RightInfo+LeftInfo);
    end
end

        
            
            
        
