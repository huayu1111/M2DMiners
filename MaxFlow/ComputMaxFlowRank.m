function [CandiRankValue,FlowValue,CutMatrix,RMatrix] = ComputMaxFlowRank(Heter_Network)
%==================================================================================
% [CandiRankValue,FlowValue,CutMatrix,RMatrix] =ComputMaxFlowRank(Heter_Network)
% Heter_Network--phenome-microRNAome network
% CandiRankValue--Candidate microRNAs ranked list
% detailed information see max-flow function
%==================================================================================
[FlowValue, CutMatrix, RMatrix,FMatrix]=max_flow(Heter_Network,size(Heter_Network,1)-1,size(Heter_Network,1));
FMatrix=full(FMatrix);
NonzeroIndex=find(FMatrix(:,539)~=0);
if(isempty(NonzeroIndex))
    CandiRankValue=[];
else
    CandiValue=FMatrix(NonzeroIndex,539);
    CandiValue=[NonzeroIndex,CandiValue];
    [~,IX]=sort(CandiValue(:,2));
    CandiRankValue=CandiValue(IX,:);
end
end

