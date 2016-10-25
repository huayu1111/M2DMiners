function [FinalCandiRankList,AUC]=MainMaxFlow()
%====================================================================================================
%== Building Phenome-MicroRNAome network and prioritizating the disease associated microRNA 
%== based on information flow maximization algorithm 
%== Data_Path: input data folder path
%== Proc_Path: program folder path
%====================================================================================================
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\MaxFlow\Data\Version2010\';
MiFunField=0.4;
DiPheField=0.1;
Alta=3;
Beta=5;
Gamma=10000;

%== building microRNA functional similarity network and return it 
MiFunSim=MiFunSimNetwork(DataPath,'mifunsim.xlsx');

%== integrate the target overlap score into microRNA functional similar network
% MiFunSim=InteTarOverlapIntoMiFun(DataPath,MiFunSim);

%== building microRNA-disease association network
mi2diNetwork=miRNA2diseaseNetwork(DataPath,'mi2di.xlsx',Alta,Beta,Gamma);

%== building disease semantic similar network
DiPheSim=DiPheSimNetwork(DataPath,'ID_data.xlsx','IC_data.xls');

%== integrate the phenotye simliarity into disease semantic similar network
DiPheSim=IntePheSim(DataPath,'OMIM.xlsx',DiPheSim);

%== building Phenome-MicroRNAome heterogenous network
Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);

%== leave-one-out cross-validation for disease prioritization
% [FinalCandiRankList,AUC]=LeaveOneOutValidationDiseases(DataPath,Heter_Network);

%== leave-one-out cross-validation for microRNA prioritization
% [FinalCandiRankList,AUC]=LeaveOneOutValidation(DataPath,Heter_Network);

% ab initio prediction
[FinalCandiRankList,AUC]=LeaveOneOutValidationDiseasesWithUnknownAssociation(DataPath,Heter_Network);
% count=1;
% for Num=90:10:140
%      [~,AUC_ablist(count)]=LeaveOneOutValidationDiseasesWithUnknownAssociation50(DataPath,Heter_Network,Num); %#ok<AGROW>
%      count=count+1;
% end
% permutation test
% [FinalCandiRankListPerm,AUC_Perm] = PermutationTest(DataPath,MiFunSim,DiPheSim,MiFunField,DiPheField);

% case studies
% Disease_Name='Breast Neoplasms';
% CandiRankListMaxFlow=PredictMaxFlow(DataPath,Heter_Network,Disease_Name);
end


