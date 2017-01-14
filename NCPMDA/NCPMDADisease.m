function [CandiRankList] = NCPMDADisease(candiDiseases,Disease_name,MiFunSim,DiPheSim,mi2diNetwork,miRNAIndex)
	CandiValueList=zeros(length(candiDiseases),1);
	[MnRow,~] = size(MiFunSim);
	[~,DnCol] = size(DiPheSim);
	[RnRow,RnCol] = size(mi2diNetwork);
	ncp_m = zeros(RnRow,RnCol)+10^-30;
	ncp_d = zeros(RnRow,RnCol)+10^-30;
	ncp =zeros(RnRow,RnCol);
    for j = 1:RnCol
        ncp_m(:,j) = (MiFunSim*mi2diNetwork(:,j))/norm(mi2diNetwork(:,j));
    end

	for i = 1:RnRow
		ncp_d(i,:) = (mi2diNetwork(i,:)*DiPheSim)/norm(mi2diNetwork(i,:));
	end
			
	for i = 1:MnRow
		for j = 1:DnCol
			ncp(i,j) = (ncp_m(i,j) + ncp_d(i,j))/(norm(MiFunSim(i,:))+norm(DiPheSim(:,j)));
		end
	end

	for i=1:length(candiDiseases)
		LogicIndex = strcmp(candiDiseases{i},Disease_name);
		matchIndex = find(LogicIndex,1);
		CandiValueList(i) = ncp(miRNAIndex,matchIndex);
	end
	CandiRankList(:,1)=candiDiseases;
	CandiRankList(:,3)=num2cell(CandiValueList);
end
