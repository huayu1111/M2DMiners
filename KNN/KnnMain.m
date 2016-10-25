function [CandiRankListKnn] = KnnMain()
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\Knn\Data\Version2010\';
load([DataPath,'mi2diNetwork.mat']);
DiseaseName='Breast Neoplasms';
Alta=4;
Beta=4;
[mi2diNetwork]=MiRNASetWeight(DataPath,mi2diNetwork,Alta,Beta); %#ok<NODEF>
% [FinalCandiRankList,AUC]=LeaveOneOutCrossValidationKnn(DataPath,mi2diNetwork);
[CandiRankListKnn]=PredictKnn(DataPath,DiseaseName,mi2diNetwork);
end