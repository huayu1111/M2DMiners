function [FinalCandiRankListRLSMDA,AUC] = RLSMDAMain()
clear all; clc;
%== setting initilize program input parameters
DataPath='E:\WorkDir\RLSMDA\Data\';
load([DataPath,'mi2diNetwork.mat']);
[FinalCandiRankListRLSMDA,AUC]=LeaveOneOutCrossValidationRLSMDA(DataPath,mi2diNetwork);
end