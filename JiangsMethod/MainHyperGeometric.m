function [CandiRankListHyper]=MainHyperGeometric()
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\HyperGeometric\Data\Version2010\';
DiseaseName='Breast Neoplasms';
BooleanMiRNANetwork=ConBooleanMiRNA(DataPath);
BooleanPheNetwork=ConBooleanPhe(DataPath);
% [FinalCandiRankList,AUC]=LeaveOneOutCrossValidationHyperGeometricAlgori(DataPath,BooleanMiRNANetwork,BooleanPheNetwork);
CandiRankListHyper=PredictHyperGeometricAlgori(DataPath,BooleanMiRNANetwork,BooleanPheNetwork,DiseaseName);
end