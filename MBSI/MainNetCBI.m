function [CandiRankListNetCBI]=MainNetCBI()
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\NetCBI\Data\Version2010\';
DiseaseName='Breast Neoplasms';
% [FinalCandiRankList,AUC] = LeaveOneOutCrossValidationNetCBI(DataPath);
CandiRankListNetCBI=PredictNetCBI(DataPath,DiseaseName);
end