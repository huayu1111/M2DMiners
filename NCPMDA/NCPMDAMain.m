function [FinalCandiRankListNCPMDA,AUC] = NCPMDAMain()
clear all; clc;
%== setting initilize program input parameters
DataPath='E:\WorkDir\NCPMDA\Data\';
load([DataPath,'mi2diNetwork.mat']);
[FinalCandiRankListNCPMDA,AUC]=LeaveOneOutCrossValidationNCPMDA(DataPath,mi2diNetwork);
end