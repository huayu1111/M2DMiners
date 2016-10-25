function [CandiRankListMBSI]=MainMBSI()
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\MBSI\Data\Version2010\';
DiseaseName='Breast Neoplasms';
% [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationMBSI(DataPath);
CandiRankListMBSI=PredictMBSI(DataPath,DiseaseName);
end