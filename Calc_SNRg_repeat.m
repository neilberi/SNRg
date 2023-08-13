% Script to run CalculateSNRg.m multiple times for OutputFile2 or OutputFile3

clear;
clear global;
clc;
close all;

OutputFile = 3; % Choose 2 or 3 to agree with CalculateSNRg.m

%% Run the program with appropriate settings

if (OutputFile == 2)

% OutputFile2
    for TSTrial = 1:10
        run CalculateSNRg.m
    end

else

% OutputFile3
    global StepSize_in; %#ok<*TLEV,*GVMIS>
    global StepSize_invec; %#ok<*TLEV,*GVMIS>
    StepSize_invec = [0.125*0.125*0.5 0.125*0.125 0.25*0.125, 0.5*0.125, 0.125, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5]*0.05/316;
    for it = 1:length(StepSize_invec)
        StepSize_in = StepSize_invec(round(it));
        run CalculateSNRg.m
    end 
    clear global;

end