function [CandiRankListRandomWalk]=MainRandomWalk()
clear all; clc;
%== setting initilize program input parameters
DataPath='G:\WorkDir\RandomWalk\Data\Version2010\';
DiseaseName='Breast Neoplasms';
Rate=0.9;
% [FinalCandiRankList,AUC]=LeaveOneOutCrossValidationRwrmdaAlgori(DataPath,Rate);
[CandiRankListRandomWalk]=PredictRandWalk(DataPath,DiseaseName,Rate);
end