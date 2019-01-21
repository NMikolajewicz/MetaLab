function dataExtraction_GUI(input)

% DATAEXTRACTION  extract data from graphical and tabular figures.
% User-guiding prompts are provided for clarity and ease of use
%
%   [EXTRACTEDDATA] = DATAEXTRACTION(DIRECTORY, DATANAME, FILE, SAVE, SAVEFORMAT)
%   extracts data sequentially from figures stored in DIRECTORY and saves
%   EXTRACTEDDATA to DATANAME Matlab structure (.mat). Optionally, users can specify to
%   SAVE contents of DATANAME to Microsoft Excel Spreadsheet FILE (.xlsx). Depending on
%   next step of meta analysis, SAVEFORMAT will need to be different. 
%
%       DATAEXTRACTION iterates through the entire contents of
%       DIRECTORY, and stores each datasets in data structure specified
%       by DATANAME. For each iteration of DATAEXTRACTION, DATANAME is
%       updated to include new dataset and CURRENTPROGRESS (inherent integer
%       variable created by DATAEXTRACTION) increases by increment of 1
%       with each successful iteration. CURRENTPROGRESS values ranges
%       between 1 to N, where N is number of figures found in
%       DIRECTORY.
%       Option to RETURN to previous or SKIP to next figure in directory
%       is presented as user prompt.
%
%   Input Arguments
%-----------------------
%
%   DIRECTORY
%           String specifying directory of folder containing all study figures
%           Study figures in DIRECTORY must be in .PNG format
%   DATANAME
%           String specifying name of structure (.MAT) to which results will be
%           saved
%   FILE
%           String that specifies the name of the Microsoft Excel Spreadsheet
%           to save extracted data to. (ex. 'myExtractions.xlsx'). Unnecessary
%           to save to spread sheet with each DATAEXTRACTION iteration, as it can be
%           computationally taxing. Recommended to SAVE to FILE only at
%           last iteration.
%   EXPORT    Logical (TRUE/FALSE). Specifies whether extracted data is saved
%           to Microsoft Excel Spreadsheet specified by FILE. Note that
%           extracted data is always be saved independently to MATLAB data structure,
%           regardless of EXPORT state.
%   EXPORTFORMAT
%           DATAEXTRACT will export EXTRACTEDDATA to FILE in one of two
%           formats:
%           'single'        - All datasets are saved in single sheet.
%                             Assumes one observation per data set. If
%                             multiple observations are present, max
%                             response is saved as single observation.
%                             Data saved in this format be used as input to METAANALYSIS.
%           'multiple'      - Allocates one Excel sheet per dataset under assumption
%                             that multiple observations are present per data set.
%                             Data saved in this format be used as input for FITMODEL.
%
%   updated 14.09.17

%% prepare for data extraction
try;

try; properties.iterate  = input.iterate; catch; properties.iterate = false; end;
try; properties.directory  = input.directory; catch; error('Directory not specified'); end
try; properties.dataName  = input.dataName; catch; error('export MATLAB structure not specified'); end
try; properties.export  = input.export{1}; catch; properties.export = false; end
try; properties.file  = input.file{1}; catch;  if properties.export; error('export excel file not specified'); end;  end;
try; properties.exportFormat  = input.exportFormat; catch; properties.exportFormat = 'single'; end

if ~strcmp(properties.directory(end), '\'); properties.directory = [properties.directory '\']; end; % ensure directory ends with \

properties

try load('currentProgress.mat'); catch; currentProgress = 1; end
try load(properties.dataName);
    DataFinal = extractedData; catch; disp ('New data structure created'); end;

loop = 1;
while loop == 1;
studyFigures = properties.directory;
addpath(studyFigures);
listing = dir([studyFigures '*.png']);  fileNames = listing(currentProgress).name; %list of figures in folder
totalStudies = length (listing); i = currentProgress;
studyFigure{i} = imread(fileNames);
figure; imshow(studyFigure{i});

extractionOptions = {'CONTINUE with current figure',...
    'RETURN to previous figure (delete last set of entries)',...
    'RETURN to previous figure (keep last set of entries)',...
    'SKIP this figure',...
    'RESET currentProgress to 1',...
    'EXPORT extracted data to Excel',...
    'EXIT data extracter'};
extractionChoice = menu ('Options', extractionOptions);

if extractionChoice == 1
    figure('NumberTitle', 'off', 'Name', 'select area of interest then double click to crop and zoom');
    subImage = imcrop(studyFigure{i});
    imshow(subImage);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf, 'NumberTitle', 'off', ...
        'Name', ' ');
    
    %% collect information about type of data presented
    
    figureTypes = {'Bivariate Graph',...
        'Univariate Graph (vertical)',...
        'Univariate Graph (horizontal)',...
        'Tabular'};
    figureType = menu ('Type of Input', figureTypes);
    
    if figureType == 1
        axesUnits = inputdlg({'x axis units (optional)', 'y axis units (optional)'}, 'Axes Units',  1);
        xaxisScales = {'extract graphically', 'custom input'};
        xaxisScale = menu ('X axis options', xaxisScales);
    elseif figureType == 2
        axesUnits = inputdlg({'vertical axis units (optional)'}, 'Axes Units',  1);
        xaxisScale = [];
    elseif figureType == 3
        axesUnits = inputdlg({'horizontal axis units  (optional)'}, 'Axes Units',  1);
        xaxisScale = [];
    elseif figureType == 4
        axesUnits = [];
        xaxisScale = [];
    end
    
    dataTypes = {'y Data only',...
        'y Data + y Error',...
        'Other'};
    dataType = menu('Specify what was reported', dataTypes);
    
    errorTypes = {'SEM',...
        'SD',...
        'unclear/not reported'};
    if dataType == 2; errorType = menu('Error Type', errorTypes); else; errorType = 3; end
    
    baseline = {'baseline absent', 'baseline present'};
    baselineCorrection = menu ('Reported Baseline', baseline);
    
    setCount = inputdlg('Enter number of datasets (required)', 'N datasets', 1);
    nDataSets = str2num(setCount{1});
    
    if isempty(nDataSets); error('Number of datasets was not declared. Cannot proceed.'); end
    
    for j = 1:nDataSets;
        labelPrompt{j} = ['Data set label # ' num2str(j) ': (required)'];
        trialPrompt{j} = ['Set' num2str(j) ' Treatment n: (optional)'];
        ctrlPrompt{j} = ['Set' num2str(j) ' Control n: (optional)'];
    end
    
    if nDataSets > 1
        repeatCalibration = menu ('Would you like to use the same axes calibration for all data-sets?', {'yes', 'no'});
    else
        repeatCalibration = 1;
    end
    
    setLabels = inputdlg(labelPrompt, 'Enter Dataset Labels ', 1);
    nTrials = inputdlg(trialPrompt, 'Enter Dataset Sample Sizes ', 1);
    
    if baselineCorrection == 2;
        nControlTrials = inputdlg(ctrlPrompt, 'Enter Dataset Sample Sizes (optional)', 1);
    else; for i = 1:nDataSets;  nControlTrials{i} = [];  end; end
    
    if figureType ~= 4
        
        %% custom x or y input
        if xaxisScale ~= 1 % if xaxis requires user input
            customXPrompt = {'enter x axis values, in order separated by commas (optional)'};
            customXTitle = 'custom x input ';
            customX = inputdlg(customXPrompt, customXTitle, 1);
            xaxisValues = cellfun(@str2double,(regexp(customX{1}, ',', 'split')));
        else; xaxisValues = []; end
        
        %% calibrate axes
        [A,B,C, ymax, xmax, origin] = axesCalibration(figureType);
        
        %% extract data
        Axis = []; Yaxis = []; Xaxis = [];
        for j = 1:nDataSets;
            % extract data points
            if baselineCorrection == 2;
                [xRaw{j},yRaw{j}, Axis, Yaxis, Xaxis] = DT1(subImage,A,B,C, xaxisValues, figureType, setLabels{j}, Axis, Yaxis, Xaxis, repeatCalibration, j);
                [x{j},y{j}, xCtrl{j}, yCtrl{j}] = DT2(subImage,A,B,C, Axis, Yaxis, Xaxis, xRaw{j}, yRaw{j}, setLabels{j});
                if ~isempty (xaxisValues); x{j} = xaxisValues; end;
            else; [xRaw{j},yRaw{j}, Axis, Yaxis, Xaxis] = DT1(subImage,A,B,C,  xaxisValues, figureType, setLabels{j}, Axis, Yaxis, Xaxis, repeatCalibration, j); x{j} = xRaw{j}; y{j} = yRaw{j};
                if ~isempty (xaxisValues); x{j} = xaxisValues; end;
                xCtrl{j} = []; yCtrl{j} = [];
            end
            
            % extract error bars
            if dataType == 2;
                basalEr = 0
                [xErrorTreat{j}, yErrorTreat{j}] = DT3(subImage,A,B,C, Axis, Yaxis, Xaxis,x{j},y{j}, setLabels{j}, basalEr);
                if baselineCorrection == 2
                    basalEr = 1;
                    [xErrorCtrl{j}, yErrorCtrl{j}] = DT3(subImage,A,B,C, Axis, Yaxis, Xaxis, xCtrl{j}, yCtrl{j}, setLabels{j}, basalEr);
                else; yErrorCtrl{j} = []; xErrorCtrl{j} = []; end
            else; yErrorTreat{j} = [];  yErrorCtrl{j} = []; xErrorTreat{j} = [];  xErrorCtrl{j} = [];
            end
        end
        
    elseif figureType == 4;
        for j = 1:nDataSets;
            tabularPrompt = {['dataset' num2str(j) ': ENTER PREDICTORS, separated by commas or leave blank if not applicable'],...
                ['dataset' num2str(j) ': ENTER RESPONSE, separated by commas or leave blank if not applicable'],...
                ['dataset' num2str(j) ': ENTER RESPONSE ERROR, separated by commas or leave blank if not applicable'],...
                ['dataset' num2str(j) ': ENTER BASELINE, separated by commas or leave blank if not applicable'],...
                ['dataset' num2str(j) ': ENTER BASELINE ERROR, separated by commas or leave blank if not applicable']};
            
            tabularTitle = 'input tabular data';
            tabularData = inputdlg(tabularPrompt, tabularTitle, 1);
            for i = 1:length(tabularData)
                tabD{i} = cellfun(@str2double,(regexp(tabularData{i}, ',', 'split')));
            end
            x{j} = tabD{1}; xCtrl{j} = tabD{1};
            y{j} = tabD{2}; yCtrl{j} = tabD{4};
            yErrorTreat{j} = tabD{3}; yErrorCtrl{j} = tabD{5};
        end
    end
    
    %% store data in structure
    for i = 1:nDataSets
        if figureType == 1; xCtrl{i} = x{i}; xUnit = axesUnits{1}; yUnit = axesUnits{2};
        elseif figureType == 2; x{i} = [];  xCtrl{i} = []; xUnit = []; yUnit = axesUnits{1};
        elseif figureType == 3; xUnit = axesUnits{1}; yUnit = [];
            y{i} = x{i}; x{i} = [];   yCtrl{i} = xCtrl{i}; xCtrl{i} = [];
            yErrorTreat{i} =  xErrorTreat{i};  xErrorTreat{i} =[];
            yErrorCtrl{i} = xErrorCtrl{i}; xErrorCtrl{i} = [];
        end
        dataExtract = struct ('Study', fileNames,...
            'FigureType', figureTypes{figureType},...
            'dataContent',  dataTypes{dataType},...
            'ErrorType', errorTypes{errorType},...
            'reportedBaseline', baseline{baselineCorrection},...
            'setLabel', setLabels{i},...
            'exposureUnits', xUnit,...
            'responseUnits', yUnit,...
            'exposure', x{i},...
            'xr', y{i},...
            'ser', yErrorTreat{i},...
            'nr', nTrials{i},...
            'xc', yCtrl{i},...
            'sec', yErrorCtrl{i},...
            'nc', nControlTrials{i},...
            'nDataSets', nDataSets);
        try; DataFinal = [DataFinal dataExtract]; catch; DataFinal = [dataExtract]; end
    end
    
    m = length(DataFinal);
    extractedData = DataFinal;
    extractedData(m)
    
    save(properties.dataName, 'extractedData');
    currentProgress = currentProgress + 1;
    save(['currentProgress.mat'],'currentProgress');
    
    display('Please check that extracted data is correct (above)');
    display('If unsatisfied, select "RETURN to previous figure (delete last set of entries)" option next time DATAEXTRACTION is run to make another extraction attempt on this figure');
    display(' ');
    display(['Data extraction from "' fileNames '" was successful']);
    display(['Data saved to "' properties.dataName '"']);
    display(' ');
    display(['Progress in current directory: ' num2str((round (10000 * currentProgress/totalStudies))/100) '%']);
    
    
    
elseif extractionChoice == 2
    % return to previous figure (delete last set of entires)
    currentProgress = currentProgress -1;
    save(['currentProgress.mat'],'currentProgress');
    extractedData = extractedData(1:(end-extractedData(end).nDataSets));
    save(properties.dataName, 'extractedData');
    display(' ');
    display(['Data extraction was not extracted for "' fileNames '"']);
    display(['Last entry in "'  properties.dataName '" was deleted']);
    display('next time DATAEXTRACTION is run, will start with previous figure in DIRECTORY');
    display(' ');
    
elseif extractionChoice == 3
    % return to previous figure (keep last entry)
    currentProgress = currentProgress -1;
    save(['currentProgress.mat'],'currentProgress');
    display(' ');
    display(['Data extraction was not extracted for "' fileNames '"']);
    display('next time DATAEXTRACTION is run, will start with previous figure in DIRECTORY');
    display(' ');
    
elseif extractionChoice ==4
    % skip current figure
    currentProgress = currentProgress + 1;
    save(['currentProgress.mat'],'currentProgress');
    display(' ');
    display(['Data extraction was not extracted for "' fileNames '"']);
    display('next time DATAEXTRACTION is run, will start with next figure in DIRECTORY');
    display(' ');
    
    
elseif extractionChoice ==5
    % rest current progress to 1
    currentProgress = 1;
    save(['currentProgress.mat'],'currentProgress');
    display(' ');
    display(['Progress has been reset to 1']);
    display(['Any existing datasets (ex. "' properties.dataName '") have not been cleared'])
    display(['Declare new dataset (can do this by entering different DATANAME as input) otherwise data extractions will be stored in pre-existing structure']);
    display(' ');
    
elseif extractionChoice ==6
    % save current structure containing exported data  
    properties.export = true;
elseif extractionChoice ==7
    properties.iterate = false;
end


if properties.export
    properties.file = [properties.file '.xlsx'];
    try;
        if length(extractedData) > 0
            fields2rmv = {'nDataSets'};
            for r = 1:length(fields2rmv); try; extractedData = rmfield(extractedData, fields2rmv{r}); catch; end; end
            switch properties.exportFormat
                case 'single'
                    extractedData = multiple2max(extractedData);
                    for i = 1:length(extractedData);
                      extractedData(i).ID = i; end                                    
                    SheetName = 'extractedData';
                    try; writetable(struct2table(extractedData), properties.file, 'Sheet', SheetName);
                    catch; writetable(struct2table(extractedData, 'AsArray', true), properties.file, 'Sheet', SheetName); end
                    display(['extractedData succesfully exported to "' properties.file '"']);
                case 'multiple'
                    D =   multipleExport(extractedData, properties);
                    for i = 1:length(D)
                        SheetName = ['extractedData ' num2str(i)];
                        try; writetable(struct2table(D(i).dataSet), properties.file, 'Sheet', SheetName);
                        catch; writetable(struct2table(D(i).dataSet, 'AsArray', true), properties.file, 'Sheet', SheetName); end
                    end
                    display(['extractedData succesfully exported to "' properties.file '"']);
            end
        end
    catch ME;
        getReport(ME)
        display('data structure was not saved');  end
end

if properties.iterate; 
    close all; 
else; loop = 0; 
end

end

close all;

catch; 

  msgbox({'Error while running extraction module',...
      '                                                                                    ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');  
end
%     clearvars('-except', 'extractedData', 'currentProgress');
end
%--------------------------------------------------------------------------
function [x,y, Axis, YScale, XScale] = DT1(Image,A,B,C, cX,  figureType, setLabel, Axis, Yaxis, Xaxis, repeatCalibration, progress)

% cX are custom input axis values (instead of ginput)
% modified to output data seleciton for subsequent error bar selection

figure, imshow(Image); set(gcf,'units','normalized','position',[0 0 1 1]);
title ('Calibrate Axis');

if repeatCalibration == 2 | progress == 1;
    display('spot 1');
    repeatCalibration
    progress
    Axis
    if figureType == 1;
        showPrompts = {'identify ymax', 'identify origin', 'identify xmax'};
        cal = [1 2 3];
    elseif figureType == 2;
        showPrompts = {'identify ymax', 'identify origin'};
        cal = [1 2];
        %      Axis(3,1) = 0;   Axis(3,2) = 0
    elseif figureType == 3;
        showPrompts = { 'identify origin', 'identify xmax'};
        cal = [2 3];
        %      Axis(1,1) = 0;   Axis(1,2) = 0
    end
    
    for i = 1:length(cal)
        title (['Calibrate Axis: ' showPrompts{i}]);
        Axis(cal(i),:) = ginput(1);
        hold on;
        plot(Axis(cal(i),1),Axis(cal(i),2),'rx','markersize',15);
    end
    if figureType == 2;  Axis(3,1) = Axis(2,1)+1;   Axis(3,2) = Axis(2,2);
    elseif figureType == 3; Axis(1,1) = Axis(2,1);   Axis(1,2) = Axis(2,2)+1; end
end

title (['Extract DATA for "' setLabel '" dataset (hit enter when complete)']);

YScale = abs(A(2)-B(2))/abs(Axis(1,2)-Axis(2,2));
XScale = abs(C(1)-B(1))/abs(Axis(3,1)-Axis(2,1));

for i=1:100
    try;  title (['Extract DATA at x = ' num2str(cX(i))]); catch; end;
    try; Points(i,:) = ginput(1); catch; break; end
    hold on;
    plot(Points(i,1),Points(i,2),'o');
end

Points(:,1) = (Points(:,1)- Axis(2,1))*XScale;
Points(:,2) = (( Axis(3,2)-Points(:,2)))*YScale;

x = Points(:,1)+B(1);
y = Points(:,2)+B(2);
end
%--------------------------------------------------------------------------

function [xTreat,yTreat, xCtrl, yCtrl] = DT2(Image,A,B,C, Axis, YScale, XScale, xData, yData , setLabel)

% modified to output data seleciton for subsequent error bar selection

figure, imshow(Image); set(gcf,'units','normalized','position',[0 0 1 1]);  hold on;

xDataScaled(:,1) = (xData / XScale )+ (Axis(2,1));
yDataScaled(:,2) = Axis(3,2) - (yData / YScale);

plot (xDataScaled, yDataScaled, 'x');
title (['Identify BASELINE for "' setLabel '" dataset']);

for i=1:100
    try; Points(i,:) = ginput(1); catch; break; end; hold on;
    plot(Points(i,1),Points(i,2),'o');
    title (['Identify BASELINE (corresponding to plotted data): ' num2str(i) '/'  num2str(length(xData))]);
    if i == length (xData); break; end
end
Points(:,1) = (Points(:,1)- Axis(2,1))*XScale;
Points(:,2) = (( Axis(3,2)-Points(:,2)))*YScale;

xCtrl = Points(:,1)+B(1); % ctrl x
yCtrl = Points(:,2)+B(2); % ctrl y

xTreat = xData; % treatment x
yTreat = yData; % treatment y
end

%--------------------------------------------------------------------------
function [x, y] = DT3(Image,A,B,C, Axis, YScale, XScale, xData, yData, setLabel, basalEr )

% studyFigure{i},A,B,C, Axis, Yaxis, Xaxis, x, y
% modified to take calibration from previoius data selection to permit
% error selection.

figure, imshow(Image); set(gcf,'units','normalized','position',[0 0 1 1]); hold on;
xDataScaled(:,1) = (xData / XScale )+ (Axis(2,1));
yDataScaled(:,2) = Axis(3,2) - (yData / YScale);
plot (xDataScaled, yDataScaled, 'x');
if basalEr
title (['Extract BASELINE ERROR for "' setLabel '" dataset (select upper error whisker)']);
else
  title (['Extract ERROR for "' setLabel '" dataset (select upper error whisker)']);  
end

for i=1:100
    try; Points(i,:) = ginput(1); catch; break; end;  hold on;
    plot(Points(i,1),Points(i,2),'o');
    if basalEr
    title (['extract BASELINE ERROR for "' setLabel '" dataset: '  num2str(i) '/'  num2str(length(xData)) '(select upper error whisker)']);
    else
    title (['extract ERROR for "' setLabel '" dataset: '  num2str(i) '/'  num2str(length(xData)) '(select upper error whisker)']);      
    end
    if i == length (xData);  break; end
end

Points(:,1) = (Points(:,1)- Axis(2,1))*XScale;
Points(:,2) = (( Axis(3,2)-Points(:,2)))*YScale;

x = abs(xData - Points(:,1)+B(1));
y = abs (yData -(Points(:,2)+B(2)));

end

function [A,B,C, ymax, xmax, origin] = axesCalibration(figureType)

if figureType == 1;
    calibrationPrompt = {'enter yMax (required)', 'enter xMax (required)'};
    calibrationTitle = 'Axes Calibration ';
    axesCalibrators = inputdlg(calibrationPrompt, calibrationTitle, 1);
    ymax = str2num(axesCalibrators{1});  xmax = str2num(axesCalibrators{2});
    
    originPrompt = {'enter x origin (required)', 'enter y origin (required)'};
    originTitle = 'Axes Origin';
    originCalibrators = inputdlg(originPrompt, originTitle, 1);
    xori = str2num(originCalibrators{1});  yori = str2num(originCalibrators{2});
    origin = [xori yori];
    
    A = [origin(1) ymax]; C = [xmax origin(2)]; B = origin;
elseif figureType == 2;
    
    calibrationPrompt = {'enter max y value (required)'};
    calibrationTitle = 'Axes Calibration';
    axesCalibrators = inputdlg(calibrationPrompt, calibrationTitle, 1);
    ymax = str2num(axesCalibrators{1});  xmax = [];
    
    originPrompt = { 'enter y origin (required)'};
    originTitle = 'Axes Origin';
    originCalibrators = inputdlg(originPrompt, originTitle, 1);
    yori = str2num(originCalibrators{1});
    origin = [yori yori];
    
    A = [origin(1) ymax]; C = [(origin(1)+1) origin(2)]; B = origin;
    
elseif figureType == 3;
    
    calibrationPrompt = {'enter max x value (required)'};
    calibrationTitle = 'Axes Calibration';
    axesCalibrators = inputdlg(calibrationPrompt, calibrationTitle, 1);
    xmax = str2num(axesCalibrators{1});  ymax = [];
    
    originPrompt = { 'enter x origin (required)'};
    originTitle = 'Axes Origin';
    originCalibrators = inputdlg(originPrompt, originTitle, 1);
    yori = str2num(originCalibrators{1});
    origin = [yori yori];
    
    A = [origin(1) (origin(2)+1)]; C = [xmax origin(2)]; B = origin;
end

end

function D =   multipleExport(extractedData, properties)
%restructure data for 'multiple' export option
for i = 1:length(extractedData)
    xr{i} = [extractedData(i).xr];
    ser{i} = [extractedData(i).ser];
    nr{i} = [extractedData(i).nr];
    xc{i} = [extractedData(i).xc];
    sec{i} = [extractedData(i).sec];
    nc{i} = [extractedData(i).nc];
    exposure{i} = [extractedData(i).exposure];
    for j = 1:length(xr{i});
        try; D(i).dataSet(j).Study =extractedData(i).Study; catch; end
        try; D(i).dataSet(j).FigureType =extractedData(i).FigureType;catch; end
        try; D(i).dataSet(j).dataContent =extractedData(i).dataContent;catch; end
        try; D(i).dataSet(j).ErrorType =extractedData(i).ErrorType;catch; end
        try; D(i).dataSet(j).reportedBaseline =extractedData(i).reportedBaseline;catch; end
        try; D(i).dataSet(j).setLabel =extractedData(i).setLabel;catch; end
        try; D(i).dataSet(j).exposureUnits = extractedData(i).exposureUnits; catch; end;
        try; D(i).dataSet(j).responseUnits = extractedData(i).responseUnits; catch; end;
        try; D(i).dataSet(j).exposure = exposure{i}(j);catch; D(i).dataSet(j).exposure = []; end
        try; D(i).dataSet(j).xr = xr{i}(j); catch; D(i).dataSet(j).xr = []; end
        try; D(i).dataSet(j).ser = ser{i}(j);catch; D(i).dataSet(j).ser =[]; end
        try; D(i).dataSet(j).nr = nr{i};catch;D(i).dataSet(j).nr =[]; end
        try; D(i).dataSet(j).xc = xc{i}(j);catch;D(i).dataSet(j).xc = []; end
        try; D(i).dataSet(j).sec = sec{i}(j);catch;D(i).dataSet(j).sec =[]; end
        try; D(i).dataSet(j).nc = nc{i};catch;D(i).dataSet(j).nc =[];  end        
        try; D(i).dataSet(j).ID = i;catch;D(i).dataSet(j).ID =[]; end     
    end
end

end

function extractedData = multiple2max(extractedData);
% convert multi-observational data sets to single observation (max value)
for i = 1:length(extractedData)
    exposure = []; xr = []; ser = []; xc = []; sec = []; iMax = []; yMax = [];
    exposure = [extractedData(i).exposure];
    xr = [extractedData(i).xr];
    ser = [extractedData(i).ser];
    xc = [extractedData(i).xc];
    sec = [extractedData(i).sec];
    
    [yMax, iMax]= max(xr);
    
    try; extractedData(i).exposure = exposure(iMax); catch; end;
    try; extractedData(i).xr = xr(iMax);catch; end;
    try; extractedData(i).ser = ser(iMax);catch; end;
    try; extractedData(i).xc = xc(iMax);catch; end;
    try; extractedData(i).sec = sec(iMax);catch; end;
end
end