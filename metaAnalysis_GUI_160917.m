%17.09.17 trouble shoot exclusion analysis
function metaAnalysis_GUI_160917(input)
clearvars('-except', 'input');

% METAANALYSIS  conduct general meta-analysis
%
%   conducts meta-analysis on DATACOLLECTION datasets and
%   saves results in FILE. Summary effect size is calculated as EFFECTSIZE
%   and transformed as DATATRANSFORMATION. Any necessary
%   backtransformations are computed using LOG2RAW method. Random effect
%   model is computed using interstudy variance estimated by TAU2ESTIMATOR
%   and confidence intervals are constructed using CIESTIMATOR. Studies
%   excluded from analysis are specified by EXCLUDESTUDIES and single study
%   subgroups are handled according to SINGLESTUDYSUBGORUSP. Study-level
%   forestplot uses weighting scheme specified by FORESTPLOTWEIGHTS, and
%   plotted data tranformation is specified by FORESTPLOTTRANSFORMATION.
%   Forestplot marker sizes are specified by MARKERSIZERANGE and optionally
%   effect sizes can be ordered if specified by SORTFOREST.
%
%   Input Arguments
%-----------------------
%
%   DATACOLLECTION
%           String specifying name of input data colection (see PREPDATA). (ex. 'myDataCollection')
%   FILE    String that specifies the name of the file to write results. (ex. 'myResults.xlsx')
%   EXPORT
%           logical (TRUE/FALSE). Specifies whether results are saved to FILE.
%   EFFECTSIZE
%           String specifying effect measure of interest.
%           EFFECTSIZE can be:
%           'absolute'      - evaluates absolute effect. If difference computed
%                             if negative control is reported (Default)
%           'standardized'  - evaluates standardized effect (Hedge's g).
%           'normalized'    - evaluates normalized effect.
%           'ratio'         - evaluated response ratio
%   DATATRANSFORMATION
%           String that specifies data transformation
%           DATATRANSFORMATION can be:
%           'raw'           - no transformation (Default)
%           'log'           - log10 transformation
%   LOG2RAW
%           String that specifies back-transformation method.
%           Available back-transformation methods:
%           'geometric'     - geometric mean
%                             result on raw scale will be approximation of
%                             data median
%           'naive'         - naive transformation (Default)
%                             Arithmetic mean on raw scale.
%                             Adhoc approximation of variance
%           'tailor'        - tailor series variance approximation
%                             same estimate as naive transformation. Variance
%                             approximation is handled differently.
%           'bias'          - bias correction method
%                             empirically derived correction factor for
%                             back-transformed estimate.
%   TAU2ESTIMATOR
%           String that specifies which tau2 estimator is used in
%           random effect model
%           TAU2ESTIMATOR can be:
%           'FE'            - Assume fixed effects model. Tau2 = 0.
%           'DL'            - DerSimonian-Laird estimator (Default)
%           'HS'            - Hunter-Schmidt estimator
%           'H'             - Hedges estimator
%           'HM'            - Hatung-Makambi estimator
%           'SJ'            - Sidik-Jonkman estimator
%   CIESTIMATOR
%           String that specifies which confidence interval estimator is
%           used to construct 95% confidence band. Options specify how
%           critical value is estimated
%           'z'             - z-distribution (Default)
%           't'             - t-distribution
%           'QA'            - quantile approximation
%   EXCLUDESTUDIES
%           Numeric array specifying which studies to exclude from overall analysis.
%           If no studies are excluded, EXCLUDESTUDIES = [].
%   SINGLESTUDYSUBGROUPS
%           logical (TRUE/FALSE). specifies whether single study subgroups are allowed (not recommended)
%   FORESTPLOTWEIGHTS
%           String specifying forest plot weighting scheme.
%           Weighting options are:
%           'FE'            - fixed effects weighting scheme
%           'RE'            - random effects weighting scheme
%   FORESTPLOTTRANSFORMATION
%           String that specifies data transformation
%           FORESTPLOTTRANSFORMATION can be:
%           'raw'           - no transformation (Default)
%           'log'           - log10 transformation
%   MARKERSIZERANGE
%           'numerical array. size(1,2). specifies min and max marker size.
%           *Adjust if forest plot markers are disproportional (Default [3,20]);
%   SORTFOREST
%           Logical (TRUE/FALSE) specifying whether forst plot effect sizes
%           are sorted in order.
%   CORRECTION
%           Logical (TRUE/FALSE) specifying whether Bonferroni correction is applied to 
%           obtain adjusted confidence intervals during subgroup analysis
%   WEIGHTS
%           String specifying the precision measure used for effect size
%           weighting. options:
%           'IV'            - inverse variance (default)
%           'IVS'           - inverse standardized variance
%           'MC'            - monte carlo resampling, pseudo raw data is
%                             randomly drawn from distribution defined by
%                             aggregate data statistics. Fixed effects not
%                             available for this option. 
%           

% Updated 01.11.17

try;
input                                                                       % summary of user inputs

% analysis properties
try; properties.effectSize = input.effectSize; catch; properties.effectSize = 'absolute'; end
try; properties.dataTransformation = input.dataTransformation; catch; properties.dataTransformation = 'raw'; end
try; properties.log2raw_method =input.log2raw_method; catch; properties.log2raw_method = 'naive'; end
try; properties.tau2estimator =input.tau2estimator; catch;  properties.tau2estimator = 'DL';end
try; properties.CIestimator = input.CIestimator; catch; properties.CIestimator = 'z';end
try; properties.export = input.export; catch;properties.export = false;  end
try; properties.excludeStudies = input.excludeStudies; catch; properties.excludeStudies = [];end
try; properties.forestPlotWeights = input.forestPlotWeights; catch; properties.forestPlotWeights = 'FE'; end
try; properties.forestPlotTransformation = input.forestPlotTransformation; catch; properties.forestPlotTransformation = 'raw'; end
try; properties.markerSizeRange = input.markerSizeRange; catch; properties.markerSizeRange = [2, 5]; end
try; properties.dataCollection = input.dataCollection; catch;  error('input data not specified'); end                                % specify dataset for meta analysis
try; properties.file = [input.file{1} '.xlsx']; catch; if properties.export; error('export file not specified'); end; end            % specify filename where results will be written
try; properties.sortForest = input.sortForest; catch; properties.sortForest = 0; end
try; properties.singleStudySub = input.singleStudySub; catch; properties.singleStudySub = false; end
try; properties.correction  = input.correction; catch; properties.correction = false; end; 
try; properties.weights = input.weights; catch; properties.weights = 'IV'; end


properties                                                                  % summary of final analysis properties

load(properties.dataCollection);                                            % load data (generated by PREPDATA module)

for k = 1:length(D)                                                         % iterate through all available data sets
    
    %% part 1: data extraction
    SheetName = D(k).description;                                           % dataset name
    data = dataExtraction(D(k).data, [properties.excludeStudies]);          % extact dataset
    assignin('base',...                                                     % assign dataset to workspace
        'extractedData',...
        dataExtraction(D(k).data,...
        [properties.excludeStudies]));
    display('part1: data extraction complete');
    
    %% part 2: compute effect size with appropriate transformations
    switch properties.effectSize
        case 'absolute'
            data =  AbsoluteDifference(data);                               % absolute difference/effect
        case 'normalized'
            data =  NormalizedDifference(data);                             % normalized difference
        case 'standardized'
            data =  hedgesG(data);                                          % standardized difference
        case 'ratio'
            data = respRatio(data);                                         % response ratio   
    end
    
    data = raw2log (data, properties, SheetName);                           % make logarithmic transformation if necessary
    display('part2: study-level effect size computation complete');    
    assignin('base', 'data', data);
    
    %% part 3: partition data into defined subgroups
    subgroup = subGroupPartition(data, D(k).covariates, properties);        % partition data into subgroups (if covariates present)
    display('part3: subgroup partitioning complete');
    
    %% part 4: conduct meta-analysis
    
    if properties.singleStudySub; minSize = 1; else; minSize = 2; end       % assign minimum subgroup size
    
    for i = 1:length(subgroup)                                              % iterate through each subgroup
        if length([subgroup(i).y]) >= (minSize)                             % analyse only if size of subgroup is sufficient
            S.subgroup = subgroup(i).subgroup;                              % specify subgroup name
            S = dataSummary(subgroup(i), S, properties);
            S = EffectSizeEstimator(S, properties);                         % summary effect size and related statistics
            try; FINAL = [FINAL S]; catch; FINAL = S; end                   % store results in final structure
        end;
    end
    
        % perform backtransformation if analyzed on log scale
        switch properties.dataTransformation                                   
        case 'log'; FINAL = log2raw(FINAL, properties); end
        display('part4: meta analysis complete');
        
        % generate study-level forest plot
        [singS] = studyLevelForestPlot (data, FINAL, properties, SheetName);
        
        %save results to spreadsheet
        FINAL = saveResults(properties, FINAL, SheetName, data);
        display(['Progress: ' num2str(100*k/(length(D))) '%'])              % print progress (% of datasets analysed)
        try; results{k} = FINAL; catch; results{1} = FINAL; end; clear('FINAL');
end

clear D i k r singleStudies  subgroup S minLen Sheet data singS;            % clear unnecessary variables
assignin('base','results', results);                                        % assign final results to workspace

catch e;
    disp(e)
    assignin('base','e',e);
      msgbox({'Error while running meta-analysis module',...
      '                                                                                    ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');  
end
end

%% FUNCTIONS

% generates study-level forest plot
function [singS] = studyLevelForestPlot(data, FINAL, properties, Sheet);

switch properties.weights
    case 'MC'
        offset = 2;
    case 'N'
        offset = 2;
    otherwise
        offset = 4;
end

switch properties.dataTransformation                                        % prep data according to desired transformation
    case 'raw'
        singS(1).study = {'Random Effects'}; singS(2).study = ' ';
        singS(1).ES = [FINAL(1).ESre ]; singS(2).ES = nan;
        singS(1).SE = [FINAL(1).SEre ]; singS(2).SE = nan;
        singS(1).xlo = singS(1).ES - (1.96*FINAL(1).SEre); singS(2).xlo = nan;
        singS(1).xhi = singS(1).ES + (1.96*FINAL(1).SEre); singS(2).xhi = nan;
        singS(1).Ni = [(FINAL(1).df+1)]; singS(2).Ni =nan;
        
        switch properties.weights 
            case 'MC'
            case 'N'
            otherwise
        singS(3).study = {'Fixed Effects'}; singS(4).study = ' ';
        singS(3).ES = [FINAL(1).ESfe ]; singS(4).ES = nan;
        singS(3).SE = [FINAL(1).SEfe ]; singS(4).SE = nan;
        singS(3).xlo = singS(3).ES - (1.96*FINAL(1).SEfe); singS(4).xlo = nan;
        singS(3).xhi = singS(3).ES + (1.96*FINAL(1).SEfe); singS(4).xhi = nan;
        singS(3).Ni = [(FINAL(1).df+1)]; singS(4).Ni = nan;
        end
        
        for i = 1:length(data);           
            singS(i+offset).study = num2str(data(i).ID);
            singS(i+offset).ES =  data(i).fESi;
            singS(i+offset).SE = data(i).fSEi;
            singS(i+offset).xlo = singS(i+offset).ES - (1.96*data(i).fSEi);
            singS(i+offset).xhi = singS(i+offset).ES + (1.96*data(i).fSEi);
            singS(i+offset).Ni = data(i).Ni;
        end
    case 'log'
        switch properties.forestPlotTransformation
            case 'raw'
                % not elegant, but it'll get the job done.
                
                %raw-scale random effect estimate
                singS(1).study = {'Random Effects'}; singS(2).study = ' ';
                singS(1).ES = [FINAL(1).ESre_raw ]; singS(2).ES = nan;
                singS(1).SE = nan; singS(2).SE = nan;
                singS(1).xlo =[FINAL(1).ciLre_raw]; singS(2).xlo = nan;
                singS(1).xhi = [FINAL(1).ciUre_raw]; singS(2).xhi = nan;
                singS(1).Ni = [(FINAL(1).df+1)]; singS(2).Ni =nan;
                
                switch properties.weights
                    case 'MC'
                    case 'N'
                    otherwise                    
                % raw-scale fixed effects estimate
                singS(3).study = {'Fixed Effects'}; singS(4).study = ' ';
                singS(3).ES = [FINAL(1).ESfe_raw ]; singS(4).ES = nan;
                singS(3).SE = nan; singS(4).SE = nan;
                singS(3).xlo = [FINAL(1).ciLfe_raw]; singS(4).xlo = nan;
                singS(3).xhi = [FINAL(1).ciUfe_raw]; singS(4).xhi = nan;
                singS(3).Ni = [(FINAL(1).df+1)]; singS(4).Ni = nan;
                end
                
                for i = 1:length(data);                                     % iterate through all data, assign to forest plot structure
                    singS(i+offset).study = data(i).ID;                          % study ID
                    singS(i+offset).ES =  data(i).ESi;                           % study-level effect size
                    singS(i+offset).SE = data(i).SEi;                            % study-level standard error
                    singS(i+offset).xlo = singS(i+offset).ES - (1.96*data(i).SEi);    % lower bound confidence interval (z-dist)
                    singS(i+offset).xhi = singS(i+offset).ES + (1.96*data(i).SEi);    % upper bound confidence interval (z-dist)
                    singS(i+offset).Ni = data(i).Ni;                             % study-level sample size
                end
            case 'log'
                
                % log-transformed Random effects estimate
                singS(1).study = {'Random Effects'}; singS(2).study = ' ';
                singS(1).ES = [FINAL(1).ESre ]; singS(2).ES = nan;
                singS(1).SE = [FINAL(1).SEre ]; singS(2).SE = nan;
                singS(1).xlo = singS(1).ES - (1.96*FINAL(1).SEre); singS(2).xlo = nan;
                singS(1).xhi = singS(1).ES + (1.96*FINAL(1).SEre); singS(2).xhi = nan;
                singS(1).Ni = [(FINAL(1).df+1)]; singS(2).Ni =nan;
                
                % log-transformed fixed effects estimate
                switch properties.weights
                    case 'MC'
                    case 'N'
                    otherwise
                singS(3).study = {'Fixed Effects'}; singS(4).study = ' ';
                singS(3).ES = [FINAL(1).ESfe ]; singS(4).ES = nan;
                singS(3).SE = [FINAL(1).SEfe ]; singS(4).SE = nan;
                singS(3).xlo = singS(3).ES - (1.96*FINAL(1).SEfe); singS(4).xlo = nan;
                singS(3).xhi = singS(3).ES + (1.96*FINAL(1).SEfe); singS(4).xhi = nan;
                singS(3).Ni = [(FINAL(1).df+1)]; singS(4).Ni = nan; 
                end
                
                for i = 1:length(data);                                     % iterate through all data, assign to forest plot structure
                    singS(i+offset).study = data(i).ID;                          % study ID
                    singS(i+offset).ES =  data(i).fESi;                          % study-level effect size
                    singS(i+offset).SE = data(i).fSEi;                           % study-level standard error
                    singS(i+offset).xlo = singS(i+offset).ES - (1.96*data(i).fSEi);   % lower bound confidence interval (z-dist)
                    singS(i+offset).xhi = singS(i+offset).ES + (1.96*data(i).fSEi);   % upper bound confidence interval (z-dist)
                    singS(i+offset).Ni = data(i).Ni;                             % study-level sample size
                end
        end
end

% sort study-level effect sizes (optional)
SE = [singS((offset+1):end).SE];                                                     % extract only study-level standard errors
ES = [singS((offset+1):end).ES];                                                     % extract only study-level effect sizes
if properties.sortForest                                                    % if sort selected
    [ES, sorted] = sort(ES, 'descend'); SE = SE(sorted);                    % sort study-level effect sizes and associated standard errors
    holdThis = singS(1:offset);
    for i = 1:length(sorted)
        holdThis = [holdThis singS(sorted(i)+offset)];
    end
    clear('singS'); singS = holdThis; clear('holdThis');
end

% % get marker weights according to specified weighting scheme
switch properties.forestPlotWeights                                         % specifies marker weighting scheme
    case  'FE'; w = 1 ./ (SE.^2);                                           % fixed effects                                        
    case 'RE'; w = 1 ./ ((SE.^2) + FINAL(1).t2);                       % random effects       
end

markerSizeRange = properties.markerSizeRange;                               % constraints on marker size range
height =  [1:(length(w))+offset];                                                % height of forest plot is length of weights and estimates

%plot forest plot
[minX, maxX, xlo, xhi, ESerrorX, ESerrorY, h] = plotSkeleton([singS.ES], [singS.SE],height, w, FINAL(1).df, 'n', [singS.xlo], [singS.xhi], markerSizeRange, [singS.study], offset, properties);

% ESerrorXr = ESerrorX; ESerrorYr = ESerrorY;

% apply confidence interval coverage for random effects
[cover1, minX1, maxX1] =  applyCoverage (ESerrorX, ESerrorY, minX, maxX, 'red');

% apply confidence interval coverage for fixed effects
switch properties.weights
    case 'MC'
    case 'N'
    otherwise
        ESerrorX{1}(1) = ESerrorX{3}(1); ESerrorX{1}(end) = ESerrorX{3}(end);
        [cover2, minX2, maxX2] =  applyCoverage (ESerrorX, ESerrorY, minX, maxX, 'blue');
end

gca(); title([Sheet ' forest plot (study-level)']);                         % name the plot
xlabel(['Effect Size (' properties.forestPlotTransformation ')']);          % label x axis
end

% evaluate coverage of meta-analysis model
function FINAL = rawModelEvaluation(FINAL);

%first check model coverage of raw data (i.e.  data from studies)
for i = 1:length(FINAL)                                                     % iterate through results
    ES = []; ciU =[]; ciL = []; ESi = [];
    ciLi = []; ciUi = []; coverage = []; cover = [];
    
    % assign data
    ES = FINAL(i).ES_raw;
    ciU = FINAL(i).ciU_raw;
    ciL =  FINAL(i).ciL_raw;
    ESi = [FINAL(i).ESi_raw];
    ciLi = [FINAL(i).ciLi_raw];
    ciUi = [FINAL(i).ciUi_raw];
    
    % compute coverage according to interval overlap
    for j = 1:length (ESi)
        cover(j) = 1;
        if ciL > ciUi(j); cover(j) = 0; end;
        if ciLi(j) > ciU; cover(j) = 0; end;
    end
    
    % store coverage results
    coverage = 100*sum(cover)/length(ESi);
    FINAL(i).studyCoverage_raw = coverage;
end

end

% backtransform log-transformed results
function FINAL = log2raw(FINAL, properties)
transformationMethod = properties.log2raw_method;                           % specify transformation metod
 
for i = 1:length(FINAL)                                                     % iterate through results                   
    zr = []; szr2 = []; n = [];  szr = [];
    
    % extract fixed effects estimates
    try;                                                                    
        zf = []; szf2 = [];  szf = [];
        zf = FINAL(i).ESfe;
        szf = FINAL(i).SEfe;
        szf2 = FINAL(i).SEfe^2;
    catch;
    end
    
     % extract random effects estimates
    zr = FINAL(i).ESre;
    szr = FINAL(i).SEre;
    szr2 = FINAL(i).SEre^2;
    
    n = FINAL(i).df + 1;                                                    % study-level sample size
    crit = FINAL(i).critValue;                                              % CI critical value
    if properties.correction
    adjCrit = FINAL(i).adjCritValue;                                        % bonferroni-correctedCI critical value
    end
    
    % conduct back-transformation according to specified method
    switch transformationMethod
        case 'tailor' % Higgins (2008)                                      % tailor series approximation method
            xr = 10^(zr+(szr2/2));                                          % z to x
            
            % standard error & confidence backtransform (random effects)
            varXr = (1/n)* ((10.^((2*zr)+(szr2))) * (szr2) * (1+((szr2)/2)));
            sxr = sqrt(varXr);
            xlor = xr - (crit * sxr);  xhir = xr + (crit * sxr);
            adj_xlor = xr - (adjCrit * sxr);  adj_xhir = xr + (adjCrit * sxr);
            
            % standard error & confidence backtransform (fixed effects)
            try;
                xf = 10^(zf+(szf2/2));
                varXf = (1/n)* ((10.^((2*zf)+(szf2))) * (szf2) * (1+((szf2)/2)));
                sxf = sqrt(varXf);
                xlof = xf - (crit * sxf);  xhif = xf + (crit * sxf);
                 adj_xlof = xf - (adjCrit * sxf);  adj_xhif = xf + (adjCrit * sxf);
            catch; end
            
            transMeth = 'tailor approximation';                             % store back-transformation method
            
        case 'bias'                                                         % bias corrected approximation, P.Rothery (1988)          
            
            % standard error approximation (random effects)
            xr = (10.^zr) +  ((1.15 * szr.^2).*(n-1 ./(n)));                % z to x
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);
            xlor = (10.^zlor) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            xhir = (10.^zhir) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            sxr = (xhir-xlor)/(2*crit);
   
            adj_zlor = zr - (adjCrit*szr); adj_zhir = zr + (adjCrit*szr);
            adj_xlor = (10.^adj_zlor) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            adj_xhir = (10.^adj_zhir) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            
            % standard error approximation (fixed effects)
            try;
                xf = (10.^zf) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                zlof = zf - (crit*szf); zhif = zf + (crit*szf);
                xlof = (10.^zlof) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                xhif = (10.^zhif) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                sxf = (xhif-xlof)/(2*crit);
                
                adj_zlof = zf - (adjCrit*szf); adj_zhif = zf + (adjCrit*szf);
                adj_xlof = (10.^adj_zlof) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                adj_xhif = (10.^adj_zhif) +  ((1.15 * szf.^2).*(n-1 ./(n)));
            catch; disp('problem with confidence intervals'); end 
            
            transMeth = 'bias correction';                                  % store back-transformation method
            
        case 'naive'                                                        % naive back-transformation
            
            % random effect backtransformation
            sxr = (10^((2*zr) + szr2))* (10^(szr2) - 1);                    % standard error (log to raw)
            xr = 10^(zr+(szr2/2));                                          % effect size (log to raw)
            
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);                 % confidence intervals (log)
            xlor = 10^(zlor+(szr2/2)); xhir = 10^(zhir+(szr2/2));           % confidence intervals (raw)
            
            if properties.correction
            adj_zlor = zr - (adjCrit*szr); adj_zhir = zr + (adjCrit*szr);   % bonferroni correct confidence intervals (log)
            adj_xlor = 10^(adj_zlor+(szr2/2));                              % bonferroni correct confidence intervals (raw)
            adj_xhir = 10^(adj_zhir+(szr2/2));
            end
            
            % fixed effect backtransformation
            try;
                sxf = (10^((2*zf) + szf2))* (10^(szf2) - 1);
                xf = 10^(zf+(szf2/2));
                zlof = zf - (crit*szf); zhif = zf + (crit*szf);
                xlof = 10^(zlof+(szf2/2));xhif = 10^(zhif+(szf2/2));
                if properties.correction
                adj_zlof = zf - (adjCrit*szf); adj_zhif = zf + (adjCrit*szf);
                adj_xlof = 10^(adj_zlof+(szf2/2));adj_xhif = 10^(adj_zhif+(szf2/2));
                end
            catch;  end
            transMeth = 'naive transformation';                             % store back-transformation method
            
        case 'geometric'                                                    % geometric mean method
            
            % random effects back-transformation
            xr = 10^(zr);                                                   % z to x                       
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);
            xlor = 10^(zlor);  xhir = 10^(zhir);
            sxr = (xhir-xlor)/(2*crit);           
            if properties.correction
            adj_zlor = zr - (adjCrit*szr); adj_zhir = zr + (adjCrit*szr);
            adj_xlor = 10^(adj_zlor);  adj_xhir = 10^(adj_zhir);
            end
            
            % fixed effects back-transformation
            try;
                xf = 10^(zf);
                zlof = zf - (crit*szf); zhif = zf + (crit*szf);
                xlof = 10^(zlof);  xhif = 10^(zhif);
                sxf = (xhif-xlof)/(2*crit);
                
                if properties.correction
                adj_zlof = zf - (adjCrit*szf); adj_zhif = zf + (adjCrit*szf);
                adj_xlof = 10^(adj_zlof);  adj_xhif = 10^(adj_zhif);
                end
            catch; end
            transMeth = 'geometric mean';                                   % store back-transformation method
    end
    
    % store back-transformed results (random effects)
    FINAL(i).ESre_raw = xr;                                                 % effect size (raw)
    FINAL(i).SEre_raw = sxr;                                                % standard error (raw)
    FINAL(i).ciLre_raw = xlor;                                              % lower bound confidence interval (raw)
    FINAL(i).ciUre_raw = xhir;                                              % upper bound confidence interval (raw)
    if properties.correction
    FINAL(i).adjCiLre_raw = adj_xlor;                                       % adjusted lower bound confidence interval (raw)
    FINAL(i).adjCiUre_raw = adj_xhir;                                       % adjusted upper bound confidence interval (raw)
    end
    
    % store back-transformed results (fixed effects)
    try;
        FINAL(i).ESfe_raw = xf;
        FINAL(i).SEfe_raw = sxf;
        FINAL(i).ciLfe_raw = xlof;
        FINAL(i).ciUfe_raw = xhif;
        if properties.correction
        FINAL(i).adjCiLfe_raw = adj_xlof;
        FINAL(i).adjCiUfe_raw = adj_xhif;
        end
    catch; end
    FINAL(i).log2rawMethod = transMeth;                                     % store transformation method
end

end

% returns heterogeneity statistics
function S = dataSummary(subgroup, S, properties)

% extract data
ES = [subgroup.y];                                                          % study-level effect sizes
SEM = [subgroup.se];                                                        % study-level standard errors
Ni = [subgroup.n];                                                          % study-level sample sizes
indicatorNum = [subgroup.indicatorNumber];                                  % study IDs

try; nSubgroups = round(mean([subgroup.nSubgroups]));                       % number of subgroups per given covariate
catch; nSubgroups = 1; end

df = length(ES) - 1;                                                        % degree of freedom
if df > 0                                                                   % if there is atleast 2 studies                         
   
    % heterogeneity statistics
    switch properties.weights
        case 'IV'
            W_fix = 1./(SEM.^2);
            ES_fix = sum(W_fix.*ES) / sum(W_fix);
        case 'IVS'
            CVi = abs(SEM./ES);  SEM = [];
            W_fix = 1./(CVi.^2); 
            SEM = abs(CVi * (sum((W_fix.*ES))/sum(W_fix)));
            W_fix = 1./(SEM.^2);
            ES_fix = sum(W_fix.*ES) / sum(W_fix);
        case 'N'
            W_fix = [Ni];
            ES_fix = sum(W_fix.*ES) / sum(W_fix);
        case 'MC'
            W_fix = 1./(SEM.^2);
            ES_fix = sum(W_fix.*ES) / sum(W_fix); 
    end

    Q = sum((W_fix).*((ES - ES_fix).^2));                                   % Q-statistic
    p_Q = 1-(round(1000*chi2cdf(Q,df))/1000);                               % p-value for Q-test
    H2 = Q / df;                                                            % higgins & thompson 2002, AKA Birge's ratio (Birge 1932).
    if Q > df;  I = 100*(H2 - 1)/H2;                                        % I2 statistic
    else; I = 0; end
    
    % confidence intrevals for h2 and I2
    if Q > df+1; SE_lnH = (0.5)*((log(Q)-log(df))/((sqrt(2*Q))-(sqrt((2*(df+1))-3))));
    else; SE_lnH = sqrt((1/(2*(df-1)))*(1-(1/(3*((df-1)^2))))); end
    H2_ciU = exp((log(H2) + (1.96 * SE_lnH)));
    H2_ciL = exp((log(H2) - (1.96 * SE_lnH)));
    I_ciU =  real(100*(H2_ciU - 1)/H2_ciU);
    I_ciL =  real(100*(H2_ciL - 1)/H2_ciL);
    if I_ciU < 0; I_ciU = 0; end; if isnan(I_ciU); I_ciU = 0; end
    if I_ciL < 0; I_ciL = 0; end; if isnan(I_ciL); I_ciL = 0; end
  
else;                                                                       % doesn't compute heterogeneity statistics for 1 study
    try; W_fix = 1./(SEM.^2);I = []; I_ciL = []; I_ciU = []; Q = []; H2 = [];p_Q = [];
    catch; I = []; I_ciL = []; I_ciU = []; Q = []; W_fix = []; H2 = [];  p_Q = []; end
end

% store statistics to structure
S.I2 = I;                                                                   % I2 heterogeneity statistic
S.I2lo = I_ciL;                                                             % lower bound I2 CI
S.I2hi = I_ciU;                                                             % upper bound I2 CI
S.H2 = H2;                                                                  % H2 heterogeneity statistic
S.df = df;                                                                  % degrees of freedom
S.Q = Q;                                                                    % Q heterogeneity statistic
S.Qp = p_Q;                                                                 % p-value for Q-test
S.xi = ES;                                                                  % study-level effect sizes for given subgroup
S.wi_fe = W_fix;                                                            % study-level weights for given subgroup
S.sei_fe = SEM;                                                             % study-level standard errors for given subgroup
S.ni = Ni;                                                                  % study-level sample sizes for given subgroup
S.id = indicatorNum;                                                        % study IDs for given subgroup
S.nSubgroups = nSubgroups;                                                  % total number of subgroups for given covariate

end

% partition data into subgroups, according to covariates
function subgroup = subGroupPartition(dataConstruct, covariates, properties);

nameHolder = dataConstruct; clear('dataConstruct');
if properties.singleStudySub; minSubSize = 1; else; minSubSize = 2; end    % specifies minimum subgroup size

% declare matrices
y = zeros(length(nameHolder),1 ); y(y == 0) = nan;                         % produce effect size array filled with nan
semy = zeros(length(nameHolder),1 ); semy(semy == 0) = nan;                % produce standard error array filled with nan
Ni = zeros(length(nameHolder),1 ); Ni(Ni == 0) = nan;                      % produce sample size array filled with nan

% assign data to matrices
for i = 1:length (nameHolder);                                             % populate arrays from above with available data
    indicatorNumber(i) = nameHolder(i).ID;
    if ~isempty(nameHolder(i).fESi); y(i) = nameHolder(i).fESi; end
    if ~isempty(nameHolder(i).fSEi);   semy(i) = nameHolder(i).fSEi; end
    if ~isempty(nameHolder(i).Ni); Ni(i) = nameHolder(i).Ni; end
end

allFields = fieldnames(nameHolder);                                        % extract all structure field names

assignin('base', 'allFields', allFields);                                  % assign structure field names to workspace 
assignin('base', 'covariates', covariates);                                % assign covariate names to workspace

% find location of covariates of interest
n = 1;
for i = 1:length(covariates)                                               % iterate through known covariates 
    for j = 1:length(allFields)                                            % iterate through datastructure field names
        if strcmp(covariates{i}, allFields{j});                            % identify indicies at which covariate is stored in structure
            keepIndex(n) = j; n = n + 1;                                
        end
    end
end

% remove unecesary fields in structure (only covariates and results remain)
try; keepIndex; catch; keepIndex = []; end                                 % check if covariates are present
for k = 1:length(allFields)
    if ~any(k == keepIndex)
        try; nameHolder = rmfield(nameHolder, allFields{k});               % remove non-covariate fields
        catch; display([  field2rmv{k} ' could not be removed from subgroup analysis because it does not exist']); end
    end
end
try; nameHolder = rmfield(nameHolder, 'ISR'); catch; end;

% group data by covariate subgroups
fnames = fieldnames(nameHolder);                                           % extract field names (corresponding to only covariates)
n = 1;
for i = 1:length(fnames)                                                   % iterate through each covariate    
    temp{i} = [nameHolder.(fnames{i})];                                    % extract covariate values
    temp{i} = temp{i}(~isnan(temp{i}));                                    % remove missing/not available data
    try;  temp{i} = str2num( temp{i}); catch; end                          % ensure covariate values are numerical
    uniq{i} = unique (temp{i});                                            % find unique covariate subgroups
    
    % for each subgroup, extract data.
    reported = [];
    startN = n;                                                            % start index for given covariate
     mSub = 0;                                                             % track number of subgroups per covariate
    for k = 1:length(uniq{i})                                              % iterature through each subgroup for given covariate
        included = [];
        count = 1;      
        try; headers{n} = [fnames{i} '_' num2str(uniq{i}(k))];             % name for covariate subgroup assigned
            y_sg = []; sd_sg = []; se_sg = []; n_sg = [];  
            yraw_sg = []; stdyraw_sg = []; semyraw_sg = []; id_sg = [];
            for m = 1:length(nameHolder)                                   %iterate through entire data structure
                
                if uniq{i}(k) == nameHolder(m).(fnames{i})                 % extract effect size and associate statistics if part of given subgroup.
                    cont = 1;                                              % flag specifyig inclusion of data in subgroup
                    
                    % ensure subgroup data is reported and valid
                    if isnan(y(m)); cont = 0; end;
                    if y(m) == 0; cont = 0; end;
                    if isempty(y(m)); cont = 0; end
                    if isnan(semy(m)); cont = 0; end; if isempty(semy(m)); cont = 0; end;
                    if isnan(Ni(m)); cont = 0; end; if isempty(Ni(m)); cont = 0; end
                    
                    % assign data to subgroup 
                    if cont == 1                                           % if all conditions are fullfilled, assign data.
                        reported = [reported m];                      
                        id_sg(count) = indicatorNumber(m);                
                        y_sg(count) = y(m);
                        se_sg(count) = semy(m);
                        n_sg(count) = Ni(m);
                        count = count + 1;
                        included(count) = m;                               % track which studies are assigned to subgroup
                    end
                end
            end
            
            % if subgroup is not empty, assign data to subgroup structure
            if ~isempty(y_sg) & length(y_sg) >= minSubSize
                subgroup(n).subgroup = headers{n};
                subgroup(n).indicatorNumber = id_sg;
                subgroup(n).y = y_sg;
                subgroup(n).se = se_sg;
                subgroup(n).n = n_sg;
                mSub = mSub+1;
                n = n + 1;
            else
                reported = setdiff(reported, included);                    % identify which studies did report covariate
            end
        catch ME;
            getReport(ME)
        end
    end
    
    count = 1;
    
    % create additional subgroup for given covariate that pools remainder
    % of studies that did not report covariate
    try; headers{n} = [fnames{i} '_nr'];                                   % characteristic not reported for this set of subgroups
        y_sg = []; sd_sg = []; se_sg = []; n_sg = [];  yraw_sg = []; 
        stdyraw_sg = []; semyraw_sg = []; id_sg = [];
        for m = 1:length(nameHolder)
            if ~ismember(m, reported)                                      % determine which study haven't been assigned to subgroup
                cont = 1;
                if isnan(y(m)); cont = 0; end;
                if y(m) == 0; cont = 0; end;
                if isempty(y(m)); cont = 0; end
                if isnan(semy(m)); cont = 0; end; if isempty(semy(m)); cont = 0; end;
                if isnan(Ni(m)); cont = 0; end; if isempty(Ni(m)); cont = 0; end
                
                if cont == 1                                               % if all conditions are fullfilled, assign data.
                    id_sg(count) = indicatorNumber(m);
                    y_sg(count) = y(m);
                    se_sg(count) = semy(m);
                    n_sg(count) = Ni(m);
                    count = count + 1;
                end
            end
        end
        
        % store not-reported subgroup to final subgroup structure
        if ~isempty(y_sg)                                               
            subgroup(n).subgroup = headers{n};
            subgroup(n).indicatorNumber = id_sg;
            subgroup(n).y = y_sg;
            subgroup(n).se = se_sg;
            subgroup(n).n = n_sg;
            mSub = mSub+1;
            n = n + 1;
        end
    catch ME;
        getReport(ME)
    end

    for nRange = startN:n
        subgroup(nRange).nSubgroups = mSub;                                % tracks number of subgroups per covariate, for later error rate correction
    end
end

% assign total set of data (with no subgroup partition) to subgroup structure
count = 1;
for m = 1:length(y);
    cont = 1;
    if isnan(y(m)); cont = 0; end;  if y(m) == 0; cont = 0; end
    if isnan(semy(m)); cont = 0; end
    if isnan(Ni(m)); cont = 0; end
    if cont == 1;
        idTotal(count) = indicatorNumber(m);
        yTotal(count) = y(m);
        seTotal(count) = semy(m);
        nTotal(count) = Ni(m);
        count = count + 1;
    end
end
total.subgroup = 'Total';
total.indicatorNumber = idTotal;
total.y = yTotal;
total.se = seTotal;
total.n = nTotal;
total.nSubgroups = 1;

try; subgroup = [total subgroup]; catch
    subgroup = [total]; end
end


function modelEval = modelEvaluation(ES, ciWidth, ESi, SEi);
ciU = ES+ciWidth;
ciUi = ESi ;
ciL = ES-ciWidth;
ciLi = ESi ;
cover = ones(1,length(ESi));
cover(ciL > ciUi) = 0;
cover(ciLi > ciU) = 0;
coverage = 100*sum(cover)/length(ESi);
modelEval.coverage = coverage;
end


function data = raw2log (data, properties, Sheet);
switch properties.dataTransformation
    case 'log';
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi) &&  data(i).ESi > 0;
                sx(i) = data(i).SEi;  Ni(i) = data(i).Ni;
                ESi(i) = data(i).ESi;
                
                sz(i) = sqrt(log10(((sx(i)^2 )/ (ESi(i)^2)) + 1));
                z(i) = log10((ESi(i)^2)/sqrt((sx(i)^2)+(ESi(i)^2)));
%                 Wz(i)= 1/ sz(i)^2;
                
                data(i).fSEi = sz(i) ;
%                 data(i).fWi = Wz(i);
                data(i).fESi = z(i);
            end
        end
        figure;
        ESi = ESi (ESi~=0);  z = z (z~=0);
        subplot(121); h1 = histogram(ESi); title([Sheet ' raw distribution']); ylabel('unweighted count'); xlabel('ESi, raw');
        subplot(122); h2 = histogram(z); title([Sheet ' log-transformed distribution']); ylabel ('unweighted count'); xlabel('ESi, log10 scale');
        set(h1, 'FaceColor' , 'k', 'EdgeColor', 'k');
        set(h2, 'FaceColor' , 'k', 'EdgeColor', 'k');
    case 'raw'
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi)
                z(i) = data(i).ESi;
                data(i).fSEi = data(i).SEi;
%                 data(i).fWi = data(i).Wi;
                data(i).fESi = data(i).ESi;
            end
        end
        figure;      
        z = z (z~=0);        
        h = histogram(z); title([Sheet ' raw distribution']); ylabel('unweighted count'); xlabel('ESi, raw');
        set(h, 'FaceColor' , 'k', 'EdgeColor', 'k');        
end
keeper = 1;
for i = 1:length(data);
    if ~isempty(data(i).fESi) && ~isnan(data(i).fESi) &&  data(i).fSEi~=0
        temp(keeper) = data(i);
        keeper = keeper + 1;
    end
end
data = []; data = temp; clear ('temp');

end

%% absolute difference
function data =  AbsoluteDifference(data)

for i = 1:length(data);
    if ~isempty(data(i).nc) &&   ~isempty(data(i).sec) &&  ~isempty(data(i).xc)
        
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        
        if ~isnan(nr(i)) && ~isnan(nc(i)); N(i) = nr(i) + nc(i); end
        if ~isnan(nr(i)) && isnan(nc(i)); N(i) = nr(i) + nr(i); end
        if  isnan(nr(i)) && ~isnan(nc(i)); N(i) = nc(i) + nc(i); end
        
        sdc(i) = sec(i) * sqrt(nc(i)); sdr(i) = ser(i) * sqrt(nr(i));
        sp(i) = sqrt((((nc(i)-1)*(sdc(i)^2))+((nr(i)-1)*(sdr(i)^2)))/(N(i)-2));
        se(i) = sqrt ((N(i)*(sp(i)^2))/(nr(i)*nc(i)));
        es(i) = xr(i) - xc(i); % absolute difference
%         w(i) = 1/(se(i)^2);
%         esw(i) = es(i)*w(i);
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
%         data(i).Wi = w(i);
%         data(i).ESWi = esw(i);
    else
        try;
            nr(i) = data(i).nr;
            ser(i) = data(i).ser;
            xr(i) = data(i).xr;
            N(i) = nr(i);
            sdr(i) = ser(i) * sqrt(nr(i));
            sp(i) = sdr(i);
            se(i) =sp(i) / sqrt(N(i));
            es(i) = xr(i); % absolute difference
%             w(i) = 1/(se(i)^2);
%             esw(i) = es(i)*w(i);
            
            data(i).Ni = N(i);
            data(i).SEi = se(i);
            data(i).ESi = es(i);
%             data(i).Wi = w(i);

        catch ME
            data(i).Ni = [];
            data(i).SEi = [];
            data(i).ESi = [];
%             data(i).Wi = [];

        end
    end
    
end
end

%% normalized difference
function data =  NormalizedDifference(data)

for i = 1:length(data);
    try;
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        
        if ~isnan(nr(i)) && ~isnan(nc(i)); N(i) = nr(i) + nc(i); end
        if ~isnan(nr(i)) && isnan(nc(i)); N(i) = nr(i) + nr(i); end
        if  isnan(nr(i)) && ~isnan(nc(i)); N(i) = nc(i) + nc(i); end
        
        sdc(i) = (sec(i) * sqrt(nc(i)))/xc(i);
        sdr(i) = (ser(i) * sqrt(nr(i)))/xr(i); %noramlized SD
        se(i) = sqrt (((sdc(i)^2)/nc(i)) + ((sdr(i)^2)/nr(i))); % normalized SE
        es(i) = (xr(i) - xc(i))/xc(i); % normalized difference
%         w(i) = 1/(se(i)^2);
%         esw(i) = es(i)*w(i);
        
        if isinf(es(i)); es(i) = nan(); end
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
%         data(i).Wi = w(i);
        
    catch ME;
        
    end
end
end

%standard difference according to Hedges
function data =  hedgesG(data)

for i = 1:length(data);
    try;
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        
        if ~isnan(nr(i)) && ~isnan(nc(i)); N(i) = nr(i) + nc(i); Np(i) = nr(i) * nc(i);end
        if ~isnan(nr(i)) && isnan(nc(i)); N(i) = nr(i) + nr(i); Np(i) = nr(i) * nr(i); end
        if  isnan(nr(i)) && ~isnan(nc(i)); N(i) = nc(i) + nc(i); Np(i) = nc(i) * nc(i);end
        
        sdc(i) = (sec(i) * sqrt(nc(i)));
        sdr(i) = (ser(i) * sqrt(nr(i)));
        sp(i) = sqrt((((nc(i)-1)*(sdc(i)^2))+((nr(i)-1)*(sdr(i)^2)))/(N(i)-2));
        es(i) = ((xr(i)-xc(i))/sp(i)) * (1-(3/((4*N(i))-9)));
        se(i) = sqrt((N(i)/Np(i))+((es(i)^2)/(2*(N(i)-3.94))));
%         
%         w(i) = 1/(se(i)^2);
%         esw(i) = es(i)*w(i);
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
%         data(i).Wi = w(i);
    catch ME;
    end
end
end



%% data extraction from structure
function [data] = dataExtraction(data, exStudies);

n = 1;
xr = []; nr = []; ser = []; xc = []; sec = []; nc = [];
for i = 1:length(data)
    if ~isempty(exStudies); exclude = any(exStudies == data(i).indicatorNumber); else; exclude = 0; end
    if ~isempty(data(i).nr) &&   ~isempty(data(i).ser) &&  ~isempty(data(i).xr) && ~isnan(data(i).nr) &&   ~isnan(data(i).ser) &&  ~isnan(data(i).xr) && exclude == 0;
        nr(n) = data(i).nr;
        ser(n) = data(i).ser;
        xr(n) = data(i).xr;
        data(i).xr = xr(n);
        data(i).nr = nr(n);
        data(i).ser = ser(n);
    end
    try;
        if ~isempty(data(i).xc) &&   ~isempty(data(i).nc) &&  ~isempty(data(i).sec) && ~isnan(data(i).xc) &&   ~isnan(data(i).nc) &&  ~isnan(data(i).sec) && exclude == 0;
            xc(n) = data(i).xc;
            try; nc(n) = data(i).nc; catch; nc(n) = nr(n); end
            sec(n) = data(i).sec;
            data(i).xc = xc(n);
            data(i).sec = sec(n);
            data(i).nc = nc(n);
        else
            data(i).xc = [];
            data(i).sec = [];
            data(i).nc = [];
        end
    catch;
        data(i).xc = [];
        data(i).sec = [];
        data(i).nc = [];
    end
    n = n+1;
end
end

function heterogeneity = tau2Estimator(S, properties)
estimator = properties.tau2estimator;

ESi = [S.xi];
n = length(ESi);
Wi = [S.wi_fe];
heterogeneity.SEi = 1./sqrt([S.wi_fe]);

vari = 1./Wi;
c = [];

switch estimator
    %% Hunter & Schmidt (HS) Estimator
    %calculates the difference between the total variance of the effect estimates and an average of the estimated within-study variances
    
    % input: n, W_fixed, ESi
    % note: n = number of studies, such that n-1 = df
    case 'HS'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        t2HS = (Q-n)/(sum(Wi));
        if Q <= n-1; t2HS = 0; end
        heterogeneity.t2 = t2HS;
        heterogeneity.t2estimator = 'Hunter_Schmidt';
        
        %% Hedges (HE) Estimator
        %calculates difference between an unweighted estimate of the total variance of the effect sizes and an unweighted estimate of the average within-study variance
        
        %input: n, ESi, vari
        %note: vari is within-study variance estimate for the ith study
    case 'H'
        ES_uw = sum(ESi) / n; %unweighted mean, i.e, arithmetic mean
        t2HE = (sum((ESi - ES_uw).^2) / (n-1)) - ((1/n) * sum(vari));
        if t2HE < 0; t2HE = 0; end % truncated at 0
        heterogeneity.t2 = t2HE;
        heterogeneity.t2estimator = 'Hedges';
        
        
        %% DerSimonian & Laird (DL) Estimator
        % moment-based method (most commonly used estimator)
        
        % input: W_fix, ESi, n
    case 'DL'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        c = sum(Wi) - ((sum((Wi).^2))/sum(Wi));
        t2DL = (Q - (n-1)) / c;
        if Q <= n-1; t2DL = 0; end % truncated at 0
        heterogeneity.t2 = t2DL;
        heterogeneity.t2estimator = 'DerSimonian_Laird';
        
        %% Hartung & Makambi (HM) Estimator
        % improvement proposed to DL estimator
        
        % input: W_fix, ESi, n
    case 'HM'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        c = sum(Wi) - ((sum((Wi).^2))/sum(Wi));
        t2HM = (Q^2)/(((2*(n-1))+Q)*c);
        if Q <= n-1; t2HM = 0; end % truncated at 0
        heterogeneity.t2 = t2HM;
        heterogeneity.t2estimator = 'Hartung_Makambi';
        
        
        %% Sidik & Jonkman (SJ) Estimator
        %based on reparametrization of the total variance in the effect size estimates
        
        %input: n, ESi
    case 'SJ'
        ES_uw = sum(ESi) / n; %unweighted mean, i.e, arithmetic mean
        t2_init = sum((ESi-ES_uw).^2)/n; %initial estimate of heterogeneity variance
        ri = vari / t2_init;
        vi = ri + 1;
        ES_v = sum((1./vi).*ESi)/sum(1./vi);
        t2SJ = (sum((1./vi).*((ESi-ES_v).^2))/(n-1));
        if t2SJ < 0; t2SJ = 0; end % truncated at 0
        heterogeneity.t2 = t2SJ;
        heterogeneity.t2estimator = 'Sidik_Jonkman';
    case 'PM'
        
        
        t2PM = 0;
        prevPM = 100;
        SEi = sqrt(vari);
        
        wCur = 1 ./ ((vari) + (t2PM^2));
        EScur = sum(wCur.*ESi)/(sum(wCur));
        
        holdThis = 1;
        while abs(t2PM - prevPM) > 0.00001
            
            prevPM(holdThis) = t2PM;
            
            wCur = 1 ./ ((vari) + (t2PM^2));
            EScur = sum(wCur.*ESi)/(sum(wCur));
            
            a = sum(wCur .*((ESi-EScur).^2));
            b = sum(wCur.*vari) - (sum((wCur.^2).*vari)/(sum(wCur)));
            c = sum(wCur) - (sum(wCur.^2)/sum(wCur));
            
            t2PM = (a-b)/c;
             if t2PM < 0; t2PM = 0; end % truncated at 0
            wCur = 1 ./ ((vari) + (t2PM^2));
            EScur = sum(wCur.*ESi)/(sum(wCur));
            holdThis = holdThis+1;
        end
%         figure; plot(1:length(prevPM), prevPM); xlabel('iteration number'); ylabel('t^2 estimate'); title('Estimation of Paule-Mandel t^2')
        heterogeneity.t2 = t2PM;
        heterogeneity.t2estimator = 'Paule-Mandel';
%         assignin('base', 'prevPM', prevPM);
%         assignin('base', 'hetergeneity', heterogeneity);
        
    case 'FE'
        heterogeneity.t2 = 0;
        heterogeneity.t2estimator = 'Fixed Effects';
end
ES_fix = (sum((Wi.*ESi))/sum(Wi));
Q = sum((Wi).*((ESi - ES_fix).^2));
heterogeneity.Q = Q;
heterogeneity.C = c;

end

function [S] = EffectSizeEstimator(S, properties)

alpha = 0.05; % intended significance cutoff

if S.nSubgroups == 1;
    mComp = 1;
else
    mComp = size(nchoosek(1:S.nSubgroups,2), 1); % max number of possible subgroup comparisons
end

if properties.correction
    if S.nSubgroups > 1; alphaCor = alpha/(mComp); % Dunn's correction (Dunn 1961) or bonferroni
    else; alphaCor = alpha;  end
else alphaCor = alpha; end


if S.df ~=0
    try heterogeneity = tau2Estimator(S, properties); catch ME
        getReport(ME)
    end
    
    xi = [S.xi];
    sei_fe = [heterogeneity.SEi];
    S.sei_fe = [];  S.sei_fe = sei_fe; 
    wi_fe = 1./(sei_fe.^2);
    k = length(xi);
    df = k-1;
    t2 = heterogeneity.t2;
    
    ES_fix = sum(xi .* wi_fe)/(sum(wi_fe));
    SE_fix = sqrt(1/sum(wi_fe));
    S.C = heterogeneity.C;
    
    try; wi_re = 1./(heterogeneity.t2 + (sei_fe.^2)); catch; wi_re = []; end
    sei_re = sqrt(1./wi_re);
    
    switch properties.weights
        case 'MC'
            [es, se] = monteCarloEstimation(xi, sei_fe, [S.ni]);
            S.ESre = es;
            S.SEre = se;
            S.ESfe = [];
            S.SEfe = [];
        case 'N'
            S.ESre = sum(xi .* [S.ni])/(sum([S.ni]));
            S.SEre = sqrt(sum(([S.ni]-1).*(sei_fe.^2))/(sum([S.ni])-length([S.ni])));
            S.ESfe = [];
            S.SEfe = [];
        otherwise
            S.ESre = sum(wi_re .* xi)/sum(wi_re);
            S.SEre = sqrt(1/sum(wi_re));
            S.ESfe = ES_fix;
            S.SEfe = SE_fix;
    end
    
    switch properties.CIestimator
        
        case 'z'
            z = norminv(1-alpha/2) ;
            S.ciWidth_re = z * sqrt(1/sum(wi_re));
            S.ciWidth_fe =  z * sqrt(1/sum(wi_fe));
            S.critValue = z;
            if properties.correction
                zAdj = norminv(1-alphaCor/2);
                S.adjCiWidth_re = zAdj * sqrt(1/sum(wi_re));
                try; S.adjCiWidth_fe =  zAdj * sqrt(1/sum(wi_fe)); catch; end
                S.adjCritValue = zAdj;
            end
            S.ciEstimator = 'zdist';
            
        case 't'
            t = tinv(1-alpha/2, df);
            S.ciWidth_re = t * sqrt(1/sum(wi_re));
            S.ciWidth_fe = t * sqrt(1/sum(wi_fe));
            S.critValue = t;
            if properties.correction
                tAdj = tinv(1-alphaCor/2, df);
                S.adjCiWidth_re = tAdj * sqrt(1/sum(wi_re));
                try; S.adjCiWidth_fe =  tAdj * sqrt(1/sum(wi_fe)); catch; end
                S.adjCritValue = tAdj;
            end
            S.ciEstimator = 'tdist';
            
        case 'QA'
            %corrected alpha not available for quantile approximation
            %method. 
            b = 2.061 + (4.902/k) + (0.756/sqrt(k)) - (0.959/log(k));
            S.ciWidth_re = b * sqrt(1/sum(wi_re));
            try; S.ciWidth_fe = b * sqrt(1/sum(wi_fe)); catch; end
            S.critValue = b;
            S.ciEstimator = 'quantileApproximation';
    end
    
    S.tauEstimator = properties.tau2estimator;
    S.t2 = t2;
    
    switch properties.weights
        case 'MC'
            modelEval_RE = modelEvaluation(S.ESre, S.ciWidth_re, xi, sei_re);
            S.coverage_re = modelEval_RE.coverage;
        case 'N'
            modelEval_RE = modelEvaluation(S.ESre, S.ciWidth_re, xi, sei_re);
            S.coverage_re = modelEval_RE.coverage;
        otherwise
            modelEval_RE = modelEvaluation(S.ESre, S.ciWidth_re, xi, sei_re);
            modelEval_FE = modelEvaluation(S.ESfe, S.ciWidth_fe, xi, sei_fe);
            
            S.coverage_re = modelEval_RE.coverage;
            S.coverage_fe = modelEval_FE.coverage;
    end

else
    S.C = [];
    switch properties.weights
        case 'MC'
            S.ESre = [S.xi];
            S.SEre = [S.sei_fe];
            S.ESfe = [];
            S.SEfe = [];
        case 'N'
            S.ESre = [S.xi];
            S.SEre = [S.sei_fe];
            S.ESfe = [];
            S.SEfe = [];
        otherwise
            S.ESre = [S.xi];
            S.SEre = [S.sei_fe];
            S.ESfe = [S.xi];
            S.SEfe = [S.sei_fe];
    end

         
    switch properties.CIestimator       
        case 'z'
            z = norminv(1-alpha/2) ;
            S.ciWidth_re = z *  S.SEre;
            S.ciWidth_fe =  z *  S.SEfe;
            S.critValue = z;
            S.ciEstimator = 'zdist';
            if properties.correction
                zAdj = norminv(1-alphaCor/2);
                S.adjCiWidth_re = zAdj * S.SEre;
                try; S.adjCiWidth_fe =  zAdj *S.SEfe; catch; end
                S.adjCritValue = zAdj;
            end
            
        case 't'
            t = tinv(1-alpha/2, df);
            S.ciWidth_re = t *  S.SEre;
            S.ciWidth_fe = t *  S.SEfe;
            S.critValue = t;
            if properties.correction
                tAdj = tinv(1-alphaCor/2, df);
                S.adjCiWidth_re = tAdj * S.SEre;
                try;  S.adjCiWidth_fe =  tAdj *  S.SEfe; catch; end
                S.adjCritValue = tAdj;
            end
            S.ciEstimator = 'tdist';
            
        case 'QA'
            b = 2.061 + (4.902/k) + (0.756/sqrt(k)) - (0.959/log(k));
            S.ciWidth_re = b *  S.SEre;
            try;  S.ciWidth_fe = b *  S.SEre; catch; end
            S.critValue = b;
            S.ciEstimator = 'quantileApproximation';
    end
    S.tauEstimator = 'not applicable';
    S.t2 = [];
    
    switch properties.weights
        case 'MC'
            S.coverage_re = [];
        case 'N'
            S.coverage_re = [];
        otherwise
            S.coverage_re = [];
            S.coverage_fe = [];
    end
end
end


function [minX, maxX, xlo, xhi, ESerrorX, ESerrorY, h] = plotSkeleton(ES, SE, height, w, df, weightingScheme, xlo, xhi, markerSizeRange, id, offset, properties)

h = figure;
dimHeight = 0.2;

wSum = sum(w);
w = w./wSum;
w = min(markerSizeRange) + (((w - min(w))./max(w)) * (max(markerSizeRange) - min(markerSizeRange)));
w(isnan(w)) = min(w);

assignin('base', 'ES', ES);
assignin('base', 'height', height);
assignin('base', 'w', w);
assignin('base', 'xlo', xlo);
assignin('base', 'id', id);

for i = 1:length(ES)
    if i == 1; color = 'r';  else; color = 'k'; end
    
    if ~isnan(xlo(i))
        xlo(i) = xlo(i); xhi(i) = xhi(i);
    else
        try;
            xlo(i) = (ES(i) - (1.96*SE(i)));
            xhi(i) =  (ES(i) + (1.96*SE(i)));
        catch; end
    end
    ESerrorX{i} = [xlo(i):(xhi(i) -xlo(i))/100 : xhi(i)];
    ESerrorY{i} = ones(1,length(ESerrorX{i})) * height(i);
    
    % random effects markers
    if i == 1;
        xdim = [xlo(i) ES(i) xhi(i)  ES(i)];
        ydim = [height(i) (height(i)-dimHeight) height(i) (height(i)+dimHeight)];
        patch(xdim, ydim, 'red'); hold on;
    end
    
    % fixed effects markers
    switch properties.weights
        case 'MC'
        case 'N'
        otherwise
            if i == 3;
                xdim = [xlo(i) ES(i) xhi(i)  ES(i)];
                ydim = [height(i) (height(i)-dimHeight) height(i) (height(i)+dimHeight)]  ;
                patch(xdim, ydim, 'blue'); hold on;
            end
    end
    
    % study level markers
    if i >offset
        plot(ES(i), height(i), 'o', 'color', color, 'MarkerFaceColor', color, 'MarkerSize', w(i-offset)); hold on;
        plot (ESerrorX{i}, ESerrorY{i},color, 'MarkerFaceColor', color);
    end
    
    
    if i == 1; minX = min(ESerrorX{i}); maxX = max(ESerrorX{i});
    else
        if minX > min(ESerrorX{i});minX = min(ESerrorX{i}); end
        if maxX < max(ESerrorX{i});maxX = max(ESerrorX{i}); end
    end
    ylim([0 max(height)+1]);
    xlabel('Effect Size');
    curAx = gca();
    set(curAx,'YTick', []);
    
end
set(gca,'XMinorTick','on');
yticks(gca, [1:i]);
yticklabels(gca, id);
ylabel('Study ID');
end

function [cover, minX, maxX] =  applyCoverage (ESerrorX, ESerrorY, minX, maxX, color)

xCover = [ESerrorX{1}(1) ESerrorX{1}(end) ESerrorX{1}(end) ESerrorX{1}(1)];
yCover = [ESerrorY{1}(1) ESerrorY{1}(1) ESerrorY{end}(1) ESerrorY{end}(1)];
cover = patch(xCover, yCover, color,...
    'EdgeColor', 'none');
set (cover, 'FaceAlpha', 0.3);
% axis equal
uistack(cover, 'bottom');

if any(minX > xCover); minX = min(xCover); end
if any(maxX < xCover); maxX = max(xCover); end

end

%% response ratio
function data = respRatio(data);
for i = 1:length(data);
    try;
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        if ~isnan(nr(i)) && ~isnan(nc(i)); N(i) = nr(i) + nc(i); end
        if ~isnan(nr(i)) && isnan(nc(i)); N(i) = nr(i) + nr(i); end
        if  isnan(nr(i)) && ~isnan(nc(i)); N(i) = nc(i) + nc(i); end
        sdc(i) = (sec(i) * sqrt(nc(i)));
        sdr(i) = (ser(i) * sqrt(nr(i)));
        es(i) = xr(i)/xc(i);
        se(i) = sqrt(((xr(i)^2)/(xc(i)^2))* (((sdr(i)^2)/(nr(i)*(xr(i)^2))) + ((sdc(i)^2)/(nc(i)*(xc(i)^2)))));
        % standard error computed using tailor approximation, under assumption cov(xr,xc) = 0.
        w(i) = 1/(se(i)^2);
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
        data(i).Wi = w(i);
    catch;
    end
end
end



function[fillhandle,msg]=jbfill_Reg(xpoints,upper,lower,color,edge,add,transparency)

%John A. Bockstege November 2006;

if nargin<7;transparency=.5;end %default is to have a transparency of .5
if nargin<6;add=1;end     %default is to add to current plot
if nargin<5;edge='k';end  %dfault edge color is black
if nargin<4;color='b';end %default color is blue

if length(upper)==length(lower) && length(lower)==length(xpoints)
    msg='';
    filled=[upper,fliplr(lower)];
    xpoints=[xpoints,fliplr(xpoints)];
    if add
        hold on
    end
    fillhandle=fill(xpoints,filled,color);%plot the data
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
    if add
        hold off
    end
else
    msg='Error: Must use the same number of points in each vector';
end
end

function [sheetName] = truncateSheet(sheetName)

sheetL = length(sheetName);
if sheetL > 31
    trunc = sheetL - 31;
else;
    trunc = 0;
end

sheetName = sheetName(1:end-trunc);
end

function [es, se] = monteCarloEstimation(z, sz, n)

% standard deviation (study-level)
sd = sz.*sqrt(n);

% monte-carlo sampling of study-level data
pseudoRaw = [];
for i = 1:length(z)
    temp = [];
            temp = z(i) + sd(i)*randn(n(i), 1);
            pseudoRaw = [pseudoRaw temp'];
end

es = mean(pseudoRaw);
sdes = std(pseudoRaw);
se = sdes/sqrt(length(z));

end



function FINAL = saveResults(properties, FINAL, SheetName, data)
        if properties.export                                                % exports results
            switch properties.weights
                case 'MC'
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id', 'ESfe', 'SEfe', 'ciWidth_fe', 'adjCritValue', 'ESfe_raw', 'SEfe_raw', 'ciLfe_raw', 'ciUfe_raw', 'adjCiLfe_raw', 'adjCiUfe_raw'};
                case 'N'
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id', 'ESfe', 'SEfe', 'ciWidth_fe', 'adjCritValue', 'ESfe_raw', 'SEfe_raw', 'ciLfe_raw', 'ciUfe_raw', 'adjCiLfe_raw', 'adjCiUfe_raw'};
                otherwise
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id'};      % remove array elements from data structure
            end
            for r = 1:length(fields2rmv);
                try; FINAL = rmfield(FINAL, fields2rmv{r}); catch; end;     % remove specified arrays
%                 try; exSens = rmfield(exSens, fields2rmv{r}); catch;  end
            end
            inputSheet = ['inputData_' SheetName{1}];                       % input sheet name
            [inputSheet] = truncateSheet(inputSheet);                       % truncate sheetname if exceeds limits
            resultsSheet = SheetName{1};                                    % results sheet name
            [resultsSheet] = truncateSheet(resultsSheet);                   % truncate sheetname if exceeds limits
            try; writetable(struct2table(data), properties.file, 'Sheet', inputSheet); catch; writetable(struct2table(data, 'AsArray', true), properties.file, 'Sheet', inputSheet); end
            try; writetable(struct2table(FINAL), properties.file, 'Sheet', resultsSheet); catch; writetable(struct2table(FINAL, 'AsArray', true), properties.file, 'Sheet', resultsSheet); end
        end
    
end


function FINAL = saveResultsv2(properties, FINAL, SheetName, data)
        if properties.export                                                % exports results
            switch properties.weights
                case 'MC'
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id', 'ESfe', 'SEfe', 'ciWidth_fe', 'adjCritValue', 'ESfe_raw', 'SEfe_raw', 'ciLfe_raw', 'ciUfe_raw', 'adjCiLfe_raw', 'adjCiUfe_raw'};
                case 'N'
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id', 'ESfe', 'SEfe', 'ciWidth_fe', 'adjCritValue', 'ESfe_raw', 'SEfe_raw', 'ciLfe_raw', 'ciUfe_raw', 'adjCiLfe_raw', 'adjCiUfe_raw'};
                otherwise
                    fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id'};      % remove array elements from data structure
            end
            for r = 1:length(fields2rmv);
                try; FINAL = rmfield(FINAL, fields2rmv{r}); catch; end;     % remove specified arrays
%                 try; exSens = rmfield(exSens, fields2rmv{r}); catch;  end
            end
            inputSheet = ['inputData_' SheetName{1}];                       % input sheet name
            [inputSheet] = truncateSheet(inputSheet);                       % truncate sheetname if exceeds limits
            resultsSheet = SheetName{1};                                    % results sheet name
            [resultsSheet] = truncateSheet(resultsSheet);                   % truncate sheetname if exceeds limits
            try; writetable(struct2table(data), properties.file, 'Sheet', inputSheet); catch; writetable(struct2table(data, 'AsArray', true), properties.file, 'Sheet', inputSheet); end
            try; writetable(struct2table(FINAL), properties.file, 'Sheet', resultsSheet); catch; writetable(struct2table(FINAL, 'AsArray', true), properties.file, 'Sheet', resultsSheet); end
        end
    
end