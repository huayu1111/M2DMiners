function [ResultsArray] = ParaOptimalWithUnknownAssociation(MiTop,DiTop,AltaTop,BetaTop)
%=============================================================================================
% [ResultsArray]=ParaOptimal(MiTop,DiTop,AltaTop,BetaTop)
%=============================================================================================
DataPath='E:\WorkDir\MaxFlow\Data\Version2010\';
ResultsArray=zeros(0,0);
counter=1;
for MiFunField=0.1:0.1:MiTop
    for DiPheField=0.1:0.1:DiTop
        for Alta=1:1:(AltaTop)
            for Beta=1:1:(BetaTop)
                load([DataPath,'Heter_Network'],'MiFunSim');
                load([DataPath,'Heter_MiDiOrigin.mat']);
                mi2diNetwork=InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,1);
                load([DataPath,'Heter_Network'],'DiPheSim');
                %== building Phenome-MicroRNAome heterogenous network
                Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);
                %== leave-one-out cross-validation
                [~,AucValue]=LeaveOneOutValidationWithUnknownAssociation(DataPath,Heter_Network);
                ResultsArray(counter,1:6)=[AucValue,MiFunField,DiPheField,Alta,Beta,1];
                counter=counter+1;
            end
        end
    end
end

for MiFunField=0.1:0.1:MiTop
    for DiPheField=0.1:0.1:DiTop
        for Alta=1:1:(AltaTop)
            for Beta=1:1:(BetaTop)
                load([DataPath,'Heter_Network'],'MiFunSim');
                load([DataPath,'Heter_MiDiOrigin.mat']);
                mi2diNetwork=10*mi2diNetwork;
                mi2diNetwork=InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,10);
                load([DataPath,'Heter_Network'],'DiPheSim');
                %== building Phenome-MicroRNAome heterogenous network
                Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);
                %== leave-one-out cross-validation
                [~,AucValue]=LeaveOneOutValidationWithUnknownAssociation(DataPath,Heter_Network);
                ResultsArray(counter,1:6)=[AucValue,MiFunField,DiPheField,Alta,Beta,10];
                counter=counter+1;
            end
        end
    end
end

for MiFunField=0.1:0.1:MiTop
    for DiPheField=0.1:0.1:DiTop
        for Alta=1:1:(AltaTop)
            for Beta=1:1:(BetaTop)
                load([DataPath,'Heter_Network'],'MiFunSim');
                load([DataPath,'Heter_MiDiOrigin.mat']);
                mi2diNetwork=100*mi2diNetwork;
                mi2diNetwork=InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,100);
                load([DataPath,'Heter_Network'],'DiPheSim');
                %== building Phenome-MicroRNAome heterogenous network
                Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);
                %== leave-one-out cross-validation
                [~,AucValue]=LeaveOneOutValidationWithUnknownAssociation(DataPath,Heter_Network);
                ResultsArray(counter,1:6)=[AucValue,MiFunField,DiPheField,Alta,Beta,100];
                counter=counter+1;
            end
        end
    end
end

for MiFunField=0.1:0.1:MiTop
    for DiPheField=0.1:0.1:DiTop
        for Alta=1:1:(AltaTop)
            for Beta=1:1:(BetaTop)
                load([DataPath,'Heter_Network'],'MiFunSim');
                load([DataPath,'Heter_MiDiOrigin.mat']);
                mi2diNetwork=1000*mi2diNetwork;
                mi2diNetwork=InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,1000);
                load([DataPath,'Heter_Network'],'DiPheSim');
                %== building Phenome-MicroRNAome heterogenous network
                Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);
                %== leave-one-out cross-validation
                [~,AucValue]=LeaveOneOutValidationWithUnknownAssociation(DataPath,Heter_Network);
                ResultsArray(counter,1:6)=[AucValue,MiFunField,DiPheField,Alta,Beta,1000];
                counter=counter+1;
            end
        end
    end
end

for MiFunField=0.1:0.1:MiTop
    for DiPheField=0.1:0.1:DiTop
        for Alta=1:1:(AltaTop)
            for Beta=1:1:(BetaTop)
                load([DataPath,'Heter_Network'],'MiFunSim');
                load([DataPath,'Heter_MiDiOrigin.mat']);
                mi2diNetwork=10000*mi2diNetwork;
                mi2diNetwork=InteMiRNASetInfo(DataPath,mi2diNetwork,Alta,Beta,10000);
                load([DataPath,'Heter_Network'],'DiPheSim');
                %== building Phenome-MicroRNAome heterogenous network
                Heter_Network=Construct_HeterNetwork(MiFunSim,DiPheSim,mi2diNetwork,MiFunField,DiPheField);
                %== leave-one-out cross-validation
                [~,AucValue]=LeaveOneOutValidationWithUnknownAssociation(DataPath,Heter_Network);
                ResultsArray(counter,1:6)=[AucValue,MiFunField,DiPheField,Alta,Beta,10000];
                counter=counter+1;
            end
        end
    end
end

