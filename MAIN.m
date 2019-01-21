%% MetaLab: Meta analysis software for basic research applications (2019)
%  N. Mikolajewicz and S.V. Komarova
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% PRESS RUN
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%%
close all; clear all; 

modules = {'Data Extraction';...
    'Fit Model';...
    'Prepare Input';...
    'Heterogeneity';...
    'Meta-analysis';...
    'Meta-regression';...
    'close MetaLab'};

figure('Name', 'MetaLab','NumberTitle','off');
try; A = imread('METALabArchitecture.png'); imshow(A); axis off; catch; end
choice = menu ('Select MetaLab Module', modules);

module = modules{choice};


switch module
    case 'Data Extraction'
        dataExtraction_UI();
    case 'Fit Model'
        fitModel_UI();
    case 'Prepare Input'
        prepData_UI();
    case 'Heterogeneity'
       h =  dataDistribution_UI();
    case 'Meta-analysis'
        metaAnalysis_UI();
    case 'Meta-regression'
        metaRegression_UI();
    case 'close MetaLab'
        close all; 
        error('MetaLab Closed');
end
