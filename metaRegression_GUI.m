function [modelStats, modelCoef, ISR_results] = metaRegression_GUI(input)

% METAREGRESSION  conduct meta-regression analysis
%
%   [RESULTS, PROPERTIES] = METAREGRESSION(DATACOLLECTION, FILE, EXPORT, MODEL,...
%   EFFECTSIZE, DATATRANSFORMATION, TAU2ESTIMATOR, ISR, ISRSEP, EXCLUDESTUDIES,...
%   MARKERSIZERANGE)conducts meta-regression analysis on
%   DATACOLLECTION datasets and EXPORTs results to FILE under assumptions of
%   fixed or random effects MODEL. Summary effect size is calculated as
%   EFFECTSIZE and transformed as DATATRANSFORMATION. Random effects model
%   is computed using interstudy variace estimated by TAU2ESTIMATOR.
%   Optionally, interstudy regression analysis is conducted for model validation
%   purposes. Studies excluded from analysis are specified by
%   EXCLUDESTUDIES and single study subgroups are handled according to
%   SINGLESTUDIES. Min and max marker sizes in final regression plot are
%   specified by MARKERSIZERANGE.
%
%   Input Arguments
%-----------------------
%
%   DATACOLLECTION
%           String specifying name of input data colection (see PREPDATA). (ex. 'myDataCollection')
%   FILE    String that specifies the name of the file to write results. (ex. 'myResults.xlsx')
%   SHEET   (Optional) String specifying name of results sheet in FILE.
%           if SHEET is not specified, default SHEET name is 'metaRegression_Results'
%   EXPORT
%           logical (TRUE/FALSE). Specifies whether results are exported to FILE.
%   MODEL   String specifying meta-analytic model. Specifically, weights are
%           calculated as inverse variance, assuming variance according to
%           'FE'            - Fixed effects model. Source of variance is
%                             sampling error (i.e. Variance = within-study variance)
%           'RE'            - Random effects model. Source of variance is
%                             sampling error and between-study heterogeneity
%                             (i.e. Variance = within-study variance + between-study variance)
%   EFFECTSIZE
%           String specifying effect measure of interest.
%           EFFECTSIZE can be:
%           'absolute'      - evaluates absolute effect. If difference computed
%                             if negative control is reported (Default)
%           'normalized'    - evaluates relative effect.
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%           'standardized'  - evaluates standardized effect.
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%           'ratio'         - evaluates response ratio
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%   DATATRANSFORMATION
%           String that specifies data transformation
%           DATATRANSFORMATION can be:
%           'raw'           - no transformation (Default)
%           'log'           - log10 transformation
%   TAU2ESTIMATOR
%           String that specifies which tau2 estimator is used in
%           random effect model. Even if fixed effects MODEL is selected,
%           random effects model is internally-computed for comparison.
%           TAU2ESTIMATOR can be:
%           'FE'            - Assume fixed effects model. Tau2 = 0.
%           'DL'            - DerSimonian-Laird estimator (Default)
%           'HS'            - Hunter-Schmidt estimator
%           'H'             - Hedges estimator
%           'HM'            - Hatung-Makambi estimator
%           'SJ'            - Sidik-Jonkman estimator
%   ISR     Logical (TRUE/FALSE) specifying whether interstudy regression
%           (ISR) is conducted. Interstudy regression fits regression-model to
%           within-study datasets and pools model parameter estimate across studies.
%           If within-study estimates are consistent with meta-regression
%           estimates of relationship between predictor and outcome, validity of
%           results are reaffirmed.
%           Interstudy regression requires ISR array as input from FILE
%           (i.e. 'ISR' header included in spreadsheet) where observations
%           drawn from the same study are specified by same value.
%           Ex.     Study               ISR
%                   'Smith 2017'        1
%                   'Smith 2017'        1
%                   'Smith 2017'        1
%                   'Foster 2017'       2
%           In this example, meta-regression model is fit to all 4
%           observations while ISR model is fit to 3 observations from
%           'Smith 2017' dataset only.
%   ISRSEP
%           Logical (TRUE/FALSE). If TRUE, ISR results are plotted separately
%           from meta-regression results. Otherwise, if FALSE, ISR results are
%           overlayed on meta-regression results.
%   EXCLUDESTUDIES
%           Numeric array specifying which studies to exclude from overall analysis,
%           according to study ID (included as header in FILE). study ID's
%           must be must be positive intergers.
%           If no studies are excluded, EXCLUDESTUDIES = [].
%   MARKERSIZERANGE
%           numerical array of size(1,2). specifies min and max marker sizes
%           for regression plot. If left empty (i.e., []), default size ranges
%           for fixed effects MODEL are [5, 50] and default size ranges for
%           random effects MODEL are [5,15]. If default sizes ranges result
%           in disproportional representation of markers, consider altering
%           range for more satisfactory illustration.
%
% Updated 16.09.17
% input

try;
    
try; properties.dataCollection = input.dataCollection; catch; error('input data not specified'); end
try; properties.model = input.model; catch; properties.model = 'FE';  end
try; properties.effectSize = input.effectSize; catch; properties.effectSize = 'absolute'; end
try; properties.dataTransformation = input.dataTransformation; catch; properties.dataTransformation = 'raw'; end
try; properties.tau2Estimator = input.tau2Estimator; catch; properties.tau2Estimator = 'DL'; end
try; properties.ISR = input.ISR; catch; properties.ISR = false; end
try; properties.ISRsep = ~input.ISRsep; catch; properties.ISRsep = true; end
try; properties.excludeStudies = input.excludeStudies; catch; properties.excludeStudies = []; end
try; properties.markerSizeRange = input.markerSizeRange; catch; properties.markerSizeRange = [];  end
try; properties.export = input.export; catch; properties.export = false; end
try; properties.file = input.file{1}; catch; if properties.export; error('export file not specified'); end; end
try; properties.Sheet = input.Sheet{1}; catch; properties.Sheet =  'metaRegression_Results'; end;

properties
properties.mechSubgroupSelect =2;

load(properties.dataCollection )

for doReg = 1:length(D)
    try;
        d = D(doReg).data;
        dataSet = D(doReg).description;
        
        parNames = D(doReg).covariates;
        parNames{end+1} = 'Study-Level Outcome';
        %     clear('D');
        
        %% model construct
        %outcome
        parLegend = {'SELECT OUTCOME'; '(continuous variable)'};
        outc = menu(parLegend, parNames);
        parNames{end} = 'fESi';
        
        %continuous predictors
        parLegend{1} = 'SELECT CONTINUOUS PREDICTORS'; parLegend{2} = '(continuous variables only)';
        parNames{end+1} = 'none of the above'; inCont = []; inCat = []; counter = 1; cat = [];
        inCont = menu(parLegend, parNames);
        if inCont ~= length(parNames); pred{counter} = parNames{inCont}; counter = counter + 1; end
        while inCont ~= length(parNames)
            inCont = menu(parLegend, parNames);
            if inCont ~= length(parNames); pred{counter} = parNames{inCont}; counter = counter + 1; end
        end
        
        %categorical predictors
        parLegend{1} = 'SELECT CATEGORICAL PREDICTORS'; parLegend{2} = '(continuous variables only)';
        inCat = menu(parLegend, parNames);
        if inCat ~= length(parNames); pred{counter} = parNames{inCat};
            cat = [cat counter]; counter = counter + 1;
        end
        while inCat ~= length(parNames)
            inCat = menu(parLegend, parNames);
            if inCat ~= length(parNames); pred{counter} = parNames{inCat};
                cat = [cat counter]; counter = counter + 1;
            end
        end
        
        %% model 1
        modelParameters.predictors = pred;
        modelParameters.outcome = parNames(outc);
        modelParameters.categorical = cat
        modelParameters.interactions = ([]);
        modelParameters.mC = [];
        
        %% subgroup selection
        
        n = 1;
        try; clear('dFinal'); catch; end
        dFinal = d;
        %     for i = 1:length(d)
        %         dFinal(n) = d(i);
        %         n = n + 1;
        %         %     end;
        %     end;
        
        Dnew(doReg).data = dFinal; Dnew(doReg).description = properties.dataCollection; clear('d'); clear('dFinal');
        data = prepData(Dnew(doReg).data, properties, modelParameters);
        % for i = 1:length(data); if data(i).MechUnits == 3; data(i).mechMagnitude = data(i).mOsmDelta; end; end
        % for i = 1:length(data); if data(i).MechUnits == 4; data(i).mechMagnitude = data(i).mOsmDelta; end; end
        % for i = 1:length(data); if data(i).MechUnits == 8; data(i).mechMagnitude = 1 / data(i).mechMagnitude ; end; end
        
        %% prep data for meta-regression
        
        modelParameters = modelConstruct(modelParameters);                                                              % model construction
        [beta, xi, yi, wi, Xo, lme, ES_sum, modelStats, ds] =  metaRegress(data, properties, modelParameters);          % meta regression
        if properties.ISR & isfield(data, 'ISR');
            [ISR, ISR_results] = intrastudyRegression(data, modelParameters, modelStats.NumCoefficients, properties);       % intrastudy regression
        else ISR = []; ISR_results = [];
            if properties.ISR & ~isfield(data, 'ISR')
                display('could not conduct intrastudy regression because ISR input not specified in data set');
            end
        end
        try;
        plotResults(modelParameters, properties, lme, beta, xi, yi, wi, ES_sum, modelStats, ISR, dataSet);                       % plot results
        catch; display('encountered problem plotting results'); end
        
        modelStats.Coefficients
        modelCoef = (modelStats.Coefficients);
        if properties.export
            inputSheet  = [dataSet{1} '_metReg_INPUT'];
            outputSheet  = [dataSet{1} '_metReg_RESULTS'];
            properties.fileOUT = [properties.file '.xlsx'];
            
            % remove variables that are not exported
            fields2rmv = modelStats.noSave;
            T = table();
            for i = 1:length(fields2rmv);  T.(fields2rmv{i}) = [modelStats.(fields2rmv{i})]; end
            modelStats = rmfield(modelStats, fields2rmv); modelStats = rmfield(modelStats, 'noSave'); modelStats = rmfield(modelStats, 'Coefficients');
            
            % export model results and statistics
            tinter.analysis = ['Meta-Regression: ' dataSet{1}];
            tintra.Coefficients = ' ';
            
            try;
                struc2xls(properties.fileOUT,tinter,'Sheet', outputSheet, 'Col', 'A', 'Row', 1, 'Orientation', 'H');
                struc2xls(properties.fileOUT,modelStats,'Sheet', outputSheet, 'Col', 'A', 'Row', 2, 'Orientation', 'H');
            catch ME; getReport(ME); end
            
            assignin('base', 'modelCoef', modelCoef);
            
            % export model coefficients (output sheet)
            modelCoefNames = modelCoef.Properties.RowNames; %get table row names
            for hob = 1:length(modelCoefNames);
                modelCoefNames{hob} = erase(modelCoefNames{hob}, ')');
                modelCoefNames{hob} = erase(modelCoefNames{hob}, '(');
                metReg.(modelCoefNames{hob}) = ''; end;
            
            try;
                struc2xls(properties.fileOUT,metReg,'Sheet', outputSheet, 'Col', 'D', 'Row', 2, 'Orientation', 'H');
                writetable(modelCoef, properties.fileOUT, 'Sheet', outputSheet, 'Range', 'E1');
            catch ME; getReport(ME)
            end
            
            %export model inputs (input sheet)
            try; writetable(T, properties.fileOUT, 'Sheet', inputSheet, 'Range', 'A1');
            catch ME; getReport(ME); end
            
            try;
                predictors = [];
                for i = 1:length(ISR_results);
                    predictors{i} = ISR_results(i).predictor;
                    slope_ISR(i) = ISR_results(i).ES; SE_ISR(i) = ISR_results(i).SE;
                    slope_meta(i) = modelCoef{predictors{i}, 'Estimate'};
                    SE_meta(i) = modelCoef{predictors{i}, 'SE'};
                    esp(i) = abs(slope_ISR(i)-slope_meta(i)); sep(i) = sqrt((SE_ISR(i)^2)+(SE_meta(i)^2));
                    tp(i) = esp(i)/sep(i);
                    F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
                    pValue_MetaVsISR(i) = quad (F, tp(i), 100)*2;
                    clear('F');
                end
                Tnew = table (predictors, slope_meta, SE_meta, slope_ISR, SE_ISR, pValue_MetaVsISR);
                struc2xls(properties.fileOUT,tintra,'Sheet', outputSheet, 'Col', 'J', 'Row', 1, 'Orientation', 'H');
                writetable(Tnew, properties.fileOUT, 'Sheet', outputSheet, 'Range', 'J2');
            catch ME; getReport(ME)
                try;
                    predictors = [];
                    coefs = modelCoef.Properties.RowNames;
                    for i = 1:height(modelCoef)-1;
                        predictors{i} = coefs{i+1};
                        slope_meta(i) = modelCoef{predictors{i}, 'Estimate'};
                        SE_meta(i) = modelCoef{predictors{i}, 'SE'};
                    end
                    Tnew = table (predictors, slope_meta, SE_meta);
                    struc2xls(properties.fileOUT,tintra,'Sheet', outputSheet, 'Col', 'J', 'Row', 1, 'Orientation', 'H');
                    writetable(Tnew, properties.fileOUT, 'Sheet', outputSheet, 'Range', 'J2');
                catch ME; getReport(ME)
                end
            end
        end
        
        
    catch ME
        getReport(ME)
        display(['error was encountered with "' dataSet{1} '" while conducting meta regression'])
       
        msgbox({'Error while running meta-regression module',...
      '                                                                                    ',...
      'Ensure inputs are correct',...
      'tip: check if there is sufficient input data for specified model'},'Error', 'Error');  
    end
end

catch
          msgbox({'Error while running meta-regression module',...
      '                                                                                    ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');  
end
end

%---------------------------------------------------------------------------
%% FUNCTIONS
function [ISR, ISR_results] = intrastudyRegression(data, modelParameters, nCoef, properties)
ISR = [];
model = properties.model;

predictors = modelParameters.predictors;
categorical = modelParameters.categorical;
outcome = modelParameters.outcome;
mC = modelParameters.mC;


% data
id = [data.ISR];
[uniqID, ia, ic] = unique(id);
n = 1;

for i = 1:length(uniqID);
    ind = []; ind = find(i == ic);
    
    if length(ind) > nCoef;
        display('happening')
        T = [];  T = table;
        ISR(n).ES = [data(ind).fESi];
        
        % specify weights, according to whether FE or RE model is used
        switch model
            case 'FE'
                ISR(n).W = [data(ind).fWi];
            case 'RE'
                ISR(n).W = [data(ind).logWi_t2];
        end
        
        % specify predictors vectors
        for k = 1:length(predictors)
            ISR(n).(predictors{k}) =  [data(ind).(predictors{k})];
            T.(predictors{k}) = [data(ind).(predictors{k})]';
        end
        
        % specify  which predictors are categorical
        for k = 1:length(categorical);
            T.(predictors{categorical(k)}) = nominal([T.(predictors{categorical(k)})]);
        end
        
        T.(outcome{1}) = [data(ind).(outcome{1})]'; % specify outcome vector
        
        ISR(n).regConstruct = T;
        
        lmeFE = fitlm(ISR(n).regConstruct, mC{1}, 'weights', [ISR(n).W]); % fit linear regression curve
        yPred = [];  yPred = predict(lmeFE); % model prediction
        ISR(n).yPred = yPred;
        
        for k = 1:length(predictors)
            exFig = figure;
            h1 = plotAdjustedResponse(lmeFE, predictors{k}); title(['adjusted: ' predictors{k}]);
            x1=get(h1,'Xdata'); y1=get(h1,'Ydata'); close(exFig);
            
            xAdj = x1{1}; xPreda = x1{2};
            yAdj = y1{1}; yPreda = y1{2};
            
            ISR(n).([predictors{k} '_xAdj']) = xAdj;
            ISR(n).([predictors{k} '_xPreda']) = xPreda;
            ISR(n).([predictors{k} '_yAdj']) = yAdj;
            ISR(n).([predictors{k} '_yPreda']) = yPreda;
            
            ISR(n).([predictors{k} '_estimate']) = lmeFE.Coefficients.Estimate(k+1);
            ISR(n).([predictors{k} '_estimateSE']) = lmeFE.Coefficients.SE(k+1);
            
            ISR(n).lme = lmeFE;
            
        end
        n = n +1;
    end
end

% calculate summary effect size for all intrastudy regression results
if ~isempty(ISR);
    for i = 1:length(predictors)
        S(i).xi =  [ISR.([predictors{i} '_estimate'])];
        S(i).wi_fe =  1./ ( [ISR.([predictors{i} '_estimateSE'])] .^2);
        
        heterogeneity = tau2Estimator(S, properties.tau2Estimator);
        t2 = heterogeneity.t2;
        
        W_re = 1 ./ ( ( [ISR.([predictors{i} '_estimateSE'])] .^2) + t2);
        
        ISR_results(i).predictor = predictors{i};
        ISR_results(i).ES = sum([S(i).xi] .* W_re)/(sum(W_re)); % summary effect size
        ISR_results(i).SE = 1/(sqrt(sum(W_re)));
        ISR_results(i).df = length([ S(i).xi ]) - 1;
        ISR_results(i).tau2 = t2;
        ISR_results(i).Q = heterogeneity.Q;
        ISR_results(i).H2 =ISR_results(i).Q  /ISR_results(i).df;
        
        if ISR_results(i).Q > ISR_results(i).df;
            ISR_results(i).I2 = 100*(ISR_results(i).H2  - 1)/ISR_results(i).H2 ;
        else; ISR_results(i).I2  = 0; end
    end
else
    ISR_results = [];
end

end




function plotRegressionConfidence(ciL, ciU, predL, predU, X, beta, xAdj, yAdj, wi, predictor, outcome, model, modelStats, properties, dataSet)
upperBound_conf =ciU; lowerBound_conf = ciL;
upperBound_pred =predU; lowerBound_pred = predL;

ES_rand_upper = modelStats.ES_rand + (1.96*modelStats.SE_rand); ES_rand_lower = modelStats.ES_rand - (1.96*modelStats.SE_rand);
ES_rand_upper = ones(size(X)) * ES_rand_upper; ES_rand_lower = ones(size(X)) * ES_rand_lower;

figure;
colorCode2 = [1 0 0];  % red
colorCode1 = [0 0 0 ]; % black

[~,~]=jbfill_Reg(X,ES_rand_lower,ES_rand_upper,colorCode1,colorCode1,0,0.1); hold on;
[~,~]=jbfill_Reg(X,lowerBound_conf,upperBound_conf,colorCode2,colorCode2,0,0.3); hold on;
[~,~]=jbfill_Reg(X,lowerBound_pred,upperBound_pred,colorCode2,colorCode2,0,0.1); hold on;

plot (X, beta(2).*X + beta(1), '--', 'Color', 'k'); hold on;

wplot = markerWeighting(wi,properties,  model);

for i = 1:length(xAdj), plot(xAdj(i), yAdj(i), 'o', 'MarkerSize', wplot(i), 'MarkerEdgeColor', 'k');  hold on; end; title('mech vs release'); ylabel('ATP release'); xlabel('mech mag'); hold on;
xlabel(predictor); ylabel(outcome)
title([dataSet ': ' predictor ' vs. ' outcome]);

xlim([min(X) max(X)]);
end


function modelParameters = modelConstruct(modelParameters)
predictors = modelParameters.predictors;
outcome = modelParameters.outcome;
interactions = modelParameters.interactions;
mC{1} = [outcome{1} '~ (' predictors{1} ')'];
for i = 2:length(predictors)
    mC{1} = [mC{1} ' + (' predictors{i} ')'];
end
for i = 1:size(interactions,1)
    mC{1} = [mC{1} ' + (' predictors{interactions(i,1)} '*' predictors{interactions(i,2)} ')'];
end
modelParameters.mC = mC;
end

function  [beta, xi, yi, wi, Xo, lme, ES_sum, modelStats, ds] = metaRegress(data, properties, modelParameters)
ds = struct2table(data);
mechSubgroupSelect = properties.mechSubgroupSelect;
model = properties.model;

mC = modelParameters.mC;
predictors = modelParameters.predictors;
outcome = modelParameters.outcome;
categorical = modelParameters.categorical;

T = table;
for i = 1:length(predictors);  T.(predictors{i}) = [data.(predictors{i})]'; end
for i = 1:length(categorical);  T.(predictors{categorical(i)}) = nominal([T.(predictors{categorical(i)})]); end
T.(outcome{1}) = [data.(outcome{1})]';
for i = 1:length(data); outY(i) = data(i).(outcome{1}); end

%fit FE model
lmeFE = fitlm(T, mC{1}, 'weights', ds.fWi);

%fi RE model
lmeRE = fitlm(T, mC{1},  'weights', ds.logWi_t2);

%Model statistics for specified model
switch model
    case 'FE';
        mW = [data.fWi]; yPred = predict(lmeFE); lme = lmeFE;
    case 'RE';
        mW = [data.logWi_t2];  yPred = predict(lmeRE); lme = lmeRE;
end

modelStats = regressionStatistics(outY, yPred, [data.fWi],lme, modelParameters, model, properties);

%diagnostic figures
figure;
subplot(231); h1 = plotResiduals(lmeFE); title('FE Residuals');  h1.FaceAlpha = .5; h1.FaceColor = 'b';
xlabel('Residuals'); ylabel('Probability');
subplot(232); plot (outY, predict(lmeFE), 'bo'); xlabel ('Observed'); ylabel ('Predicted'); title('FE: Observed vs Predicted'); hold on;
plotEquality(outY, predict(lmeFE)); legend ('obs vs pred', 'line of equality', 'location', 'best');
subplot(233); h3 = plotResiduals(lmeFE,'probability', 'Color', 'b'); legend ('line of normality', 'res vs norm', 'location', 'best');
subplot(234); h4 = plotResiduals(lmeRE); title('RE residuals'); h4.FaceAlpha = .5; h4.FaceColor = 'r';
xlabel('Residuals'); ylabel('Probability');
subplot(235); plot (outY, predict(lmeRE), 'ro'); xlabel ('Observed'); ylabel ('Predicted'); title('RE: Observed vs Predicted'); hold on;
plotEquality(outY, predict(lmeRE)); legend ('obs vs pred', 'line of equality', 'location', 'best');
subplot(236); h6 = plotResiduals(lmeRE,'probability',  'Color', 'r'); legend('line of normality', 'res vs norm', 'location', 'best');

yResp = outY;
xi = [modelStats.(predictors{1})]; yi = [modelStats.(outcome{1})]; wi = [modelStats.weights];

wplot = markerWeighting(wi, properties,model);

ESi = [modelStats.(outcome{1})]; Wi = [modelStats.weights];
ES_sum = sum(ESi .* Wi)/(sum(Wi));

Xo = []; beta = [];

end


function wplot = markerWeighting(wi, properties, model)

if isempty([properties.markerSizeRange])
    switch model
        case 'FE'
            markerSizeRange = [5, 50];
        case 'RE'
            markerSizeRange = [5, 15];
    end
else
    markerSizeRange =  [properties.markerSizeRange];
end

wi = wi./(sum(wi));
wplot = min(markerSizeRange) + (((wi - min(wi))./max(wi)) * (max(markerSizeRange) - min(markerSizeRange)));
wplot(isnan(wplot)) = min(wplot);
wplot(wplot<0) = min(wplot(wplot>0));
end

function  modelStat = regressionStatistics(obs, pred, wFix, lme, modelParameters, model, properties)

predictors = modelParameters.predictors;
outcome = modelParameters.outcome;
mC = modelParameters.mC;

%% calculate residual Q score
n = 1;
for i = 1:length(obs)
    q(i) = wFix(i) * ((obs(i)-pred(i))^2);
    if ~isnan(q(i))
        qq(n) = q(i);
        n = n+1;
    end
end
Qres = sum(qq);

%% calculate heterogeneity
for i = 1:length(obs); S(i).xi = obs(i); S(i).wi_fe = wFix(i); end
heterogeneityMetaReg = tau2Estimator(S, properties.tau2Estimator);

for i = 1:length(obs); S_rand(i).xi = obs(i); S_rand(i).wi_fe = wFix(i); end
heterogeneityRandom = tau2Estimator(S_rand, properties.tau2Estimator);

%% store statistics in structure
switch model
    case 'FE'
        modelStat.model = 'Fixed Effects';
    case 'RE'
        modelStat.model = 'Random Effects';
end

modelStat.response = lme.Formula.ResponseName;
modelStat.predictors = lme.Formula.LinearPredictor;

modelStat.R2ord = lme.Rsquared.ordinary;
modelStat.R2adj =  lme.Rsquared.adjusted;

modelStat.Qtotal = heterogeneityMetaReg.Q;  %same as presented in initial meta analysis (under fixed effects assumption)
modelStat.Qtotal_df = lme.DFE;
modelStat.Qtotal_pValue = 1-abs(chi2cdf(modelStat.Qtotal, modelStat.Qtotal_df ));

if modelStat.Qtotal_pValue < 0.05
    modelStat.Qtotal_Results = 'between-study variance > 0';
else
    modelStat.Qtotal_Results = 'between-study variance = 0';
end

modelStat.Qmodel = modelStat.Qtotal - Qres;  %if sig, relationship between outcome and covariates is greater then expected due to chance
modelStat.Qmodel_df = lme.NumCoefficients - 1;
modelStat.Qmodel_pValue = 1-abs(chi2cdf(modelStat.Qmodel, modelStat.Qmodel_df ));

if modelStat.Qmodel_pValue < 0.05
    modelStat.Qmodel_Results = 'Variance explained by model > 0';
else
    modelStat.Qmodel_Results = 'Variance explained by model = 0';
end

modelStat.Qres = Qres;  %if sig, between study variance remains underexplained
modelStat.Qres_df = modelStat.Qtotal_df - modelStat.Qmodel_df;
modelStat.Qres_pValue = 1-abs(chi2cdf(modelStat.Qres, modelStat.Qres_df ));

if modelStat.Qres_pValue < 0.05;
    modelStat.Qres_Results = 'Data not consistent with model assumptions';
else
    modelStat.Qres_Results = 'Data consistent with model assumptions';
end

modelStat.Ctotal = heterogeneityRandom.C;
modelStat.t2total = heterogeneityRandom.t2;

t2res = (modelStat.Qres - modelStat.Qres_df)/modelStat.Ctotal;
if t2res < 0; t2res = 0; end
modelStat.t2res = t2res;
modelStat.Rexplained = 100 * (1-(modelStat.t2res / modelStat.t2total));
modelStat.n = lme.NumObservations;
modelStat.NumCoefficients = lme.NumCoefficients;

df =lme.DFE;
H2 =Qres  / df; %higgins & thompson 2002, AKA Birge's ratio (Birge 1932).
if Qres > df;  I = 100*(H2 - 1)/H2;
else; I = 0; end

modelStat.I2 = I;
modelStat.LogLikelihood = lme.LogLikelihood;
modelStat.SSE =   lme.SSE;
modelStat.SST =   lme.SST;
modelStat.SSR =   lme.SSR;
modelStat.MSE =   lme.MSE;


modelStat.(outcome{1}) = lme.Variables.(outcome{1});
noSave{1} = outcome{1};
noSave{2} = 'weights';
for i = 1:length(predictors)
    modelStat.(predictors{i}) = lme.Variables.(predictors{i});
    noSave{length(noSave)+1} = predictors{i};
end
modelStat.noSave = noSave;
modelStat.weights = lme.ObservationInfo.Weights;

ESi = [modelStat.(outcome{1})]';
sei_fix = sqrt(1./wFix);
Wi_r =  1./(heterogeneityRandom.t2 + (sei_fix.^2));
modelStat.ES_fix =  sum(ESi .* wFix)/(sum(wFix));
modelStat.SE_fix =  sqrt(1/sum(wFix));
modelStat.ES_rand = sum(ESi .* Wi_r)/(sum(Wi_r));
modelStat.SE_rand = sqrt(1/sum(Wi_r));

modelStat.Coefficients = lme.Coefficients;

end

function datav2 = prepData(testSet, properties, modelParameters)

excludeStudies = properties.excludeStudies;
effectSize = properties.effectSize;

predictors = modelParameters.predictors;
outcome = modelParameters.outcome;

Sheet= []; FINAL_uniq = [];

FINAL_uniq = testSet;
[data] = dataExtraction(FINAL_uniq, excludeStudies);

switch properties.effectSize
    case 'absolute'
        data =  AbsoluteDifference(data);                               % absolute difference
    case 'normalized'
        data =  NormalizedDifference(data);                             % normalized difference
    case 'standardized'
        data =  hedgesG(data);                                            % standardized difference
    case 'ratio'
        data = respRatio(data);
end

data = raw2log (data, properties);
S.xi = [data.fESi]; S.wi_fe = [data.fWi];
heterogeneity = tau2Estimator(S, properties.tau2Estimator);

for i = 1:length(data);
    data(i).tau = heterogeneity.t2;
    data(i).Q =  heterogeneity.Q;
    
    if isnumeric(data(i).fWi) && ~isempty(data(i).fWi) && isnumeric(heterogeneity.t2) && ~isempty(heterogeneity.t2);
        data(i).logWi_t2 = 1 / ((1 / (data(i).fWi)) + heterogeneity.t2);
    else
        data(i).logWi_t2 = [];
    end
end

%remove empty, NaN fields for variables of interest.
for i = 1:length(data);
    incomplete(i) = 0;
    for j = 1:length(predictors)
        if  isempty(data(i).(predictors{j})) | ...
                isnan(data(i).(predictors{j})) | ...
                ~isnumeric(data(i).(predictors{j}))
            incomplete(i) = 1;
        end
    end
    if  isempty(data(i).fWi) | ...
            isnan(data(i).fWi) | ...
            ~isnumeric(data(i).fWi)
        incomplete(i) = 1;
    end
    if  isempty(data(i).logWi_t2) | ...
            isnan(data(i).logWi_t2) | ...
            ~isnumeric(data(i).logWi_t2)
        incomplete(i) = 1;
        
    end
    if  isempty(data(i).(outcome{1})) | ...
            isnan(data(i).(outcome{1})) | ...
            ~isnumeric(data(i).(outcome{1}))
        incomplete(i) = 1;
    end
end

datav2 = data(~incomplete);

if length(datav2)<3
    error('Not enough input observations to conduct meta-regression analysis')
end

end

function data = raw2log (data, properties);
switch properties.dataTransformation
    case 'log'
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi) &&  data(i).ESi > 0;
                sx(i) = data(i).SEi;  Ni(i) = data(i).Ni;
                ESi(i) = data(i).ESi;
                
                sz(i) = sqrt(log10(((sx(i)^2 )/ (ESi(i)^2)) + 1));
                z(i) = log10((ESi(i)^2)/sqrt((sx(i)^2)+(ESi(i)^2)));
                Wz(i)= 1/ sz(i)^2;
                
                data(i).fSEi = sz(i) ;
                data(i).fWi = Wz(i);
                data(i).fESi = z(i);
            end
        end
    case 'raw'
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi)
                z(i) = data(i).ESi;
                data(i).fSEi = data(i).SEi;
                data(i).fWi = data(i).Wi;
                data(i).fESi = data(i).ESi;
            end
        end
end
end


function heterogeneity = tau2Estimator(S, estimator)

ESi = [S.xi];
Wi = [S.wi_fe];
[ESi, Wi] = prepareCurveData(ESi, Wi);
n = length(ESi);
vari = 1./Wi;
c = [];

switch estimator
    %% Hunter & Schmidt (HS) Estimator
    %calculates the difference between the total variance of the effect estimates and an average of the estimated within-study variances
    case 'HS'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        t2HS = (Q-n)/(sum(Wi));
        if Q <= n-1; t2HS = 0; end
        heterogeneity.t2 = t2HS;
        heterogeneiry.t2estimator = 'Hunter_Schmidt';
        
        %% Hedges (HE) Estimator
        %calculates difference between an unweighted estimate of the total variance of the effect sizes and an unweighted estimate of the average within-study variance
    case 'H'
        ES_uw = sum(ESi) / n; %unweighted mean, i.e, arithmetic mean
        t2HE = (sum((ESi - ES_uw).^2) / (n-1)) - ((1/n) * sum(vari));
        if t2HE < 0; t2HE = 0; end % truncated at 0
        heterogeneity.t2 = t2HE;
        heterogeneiry.t2estimator = 'Hedges';
        
        
        %% DerSimonian & Laird (DL) Estimator
        % moment-based method (most commonly used estimator)
    case 'DL'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        c = sum(Wi) - ((sum((Wi).^2))/sum(Wi));
        t2DL = (Q - (n-1)) / c;
        if Q <= n-1; t2DL = 0; end % truncated at 0
        heterogeneity.t2 = t2DL;
        heterogeneiry.t2estimator = 'DerSimonian_Laird';
        
        %% Hartung & Makambi (HM) Estimator
    case 'HM'
        ES_fix = (sum((Wi.*ESi))/sum(Wi));
        Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
        c = sum(Wi) - ((sum((Wi).^2))/sum(Wi));
        t2HM = (Q^2)/(((2*(n-1))+Q)*c);
        if Q <= n-1; t2HM = 0; end % truncated at 0
        heterogeneity.t2 = t2HM;
        heterogeneiry.t2estimator = 'Hartung_Makambi';
        
        
        %% Sidik & Jonkman (SJ) Estimator
    case 'SJ'
        ES_uw = sum(ESi) / n; %unweighted mean, i.e, arithmetic mean
        t2_init = sum((ESi-ES_uw).^2)/n; %initial estimate of heterogeneity variance
        ri = vari / t2_init;
        vi = ri + 1;
        ES_v = sum((1./vi).*ESi)/sum(1./vi);
        t2SJ = (sum((1./vi).*((ESi-ES_v).^2))/(n-1));
        if t2SJ < 0; t2SJ = 0; end % truncated at 0
        heterogeneity.t2 = t2SJ;
        heterogeneiry.t2estimator = 'Sidik_Jonkman';
        
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
        heterogeneiry.t2estimator = 'Fixed Effects';
end

heterogeneity.C = c;
ES_fix = (sum((Wi.*ESi))/sum(Wi));
Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic
heterogeneity.Q = Q;

end


%% data extraction from structure
function [data] = dataExtraction(data, exStudies);
n = 1;
xr = []; nr = []; ser = []; xc = []; sec = []; nc = [];
for i = 1:length(data)
    if ~isempty(exStudies); exclude = any(exStudies == data(i).ISR); else; exclude = 0; end
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
        es(i) = xr(i) - xc(i);
        w(i) = 1/(se(i)^2);
        esw(i) = es(i)*w(i);
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
        data(i).Wi = w(i);
        data(i).ESWi = esw(i);
    else
        try;
            nr(i) = data(i).nr;
            ser(i) = data(i).ser;
            xr(i) = data(i).xr;
            N(i) = nr(i);
            sdr(i) = ser(i) * sqrt(nr(i));
            sp(i) = sdr(i);
            se(i) =sp(i) / sqrt(N(i));
            es(i) = xr(i);
            w(i) = 1/(se(i)^2);
            esw(i) = es(i)*w(i);
            
            data(i).Ni = N(i);
            data(i).SEi = se(i);
            data(i).ESi = es(i);
            data(i).Wi = w(i);
        catch ME
            data(i).Ni = [];
            data(i).SEi = [];
            data(i).ESi = [];
            data(i).Wi = [];
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
        w(i) = 1/(se(i)^2);
        esw(i) = es(i)*w(i);
        
        if isinf(es(i)); es(i) = nan(); end
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
        data(i).Wi = w(i);
        
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
        
        w(i) = 1/(se(i)^2);
        esw(i) = es(i)*w(i);
        
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
        data(i).Wi = w(i);
    catch ME;
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

function  plotEquality(obs, pred);

oMax = max(obs); pMax = max(pred);
oMin = min(obs); pMin = min(pred);

if oMax>pMax; ma = oMax; else; ma = pMax; end;
if oMin<pMin; mi = oMin; else; mi = pMin; end;

x = linspace(mi, ma);

plot (x,x, 'k--');
end

function plotResults (modelParameters, properties, lme, beta, xi, yi, wi, ES_sum, modelStats, ISR, dataSet)
% plotResults (modelParameters, properties, lme, beta, xi, yi, wi, Xo, ES_sum, modelStats, ds, ISR, ISR_results)

predictors = modelParameters.predictors; outcome = modelParameters.outcome;
categorical = modelParameters.categorical; interactions = modelParameters.interactions;

if length(predictors) ==1
    exFig = figure;
    h1 = plotAdjustedResponse(lme, predictors{1}); title(['adjusted: ' predictors{1}]);
    x1=get(h1,'Xdata'); y1=get(h1,'Ydata'); close(exFig);
    
    xAdj = x1{1}; xPreda = x1{2}; yAdj = y1{1}; yPreda = y1{2};
    predCurve = fitlm(xPreda, yPreda);
    beta(1) = predCurve.Coefficients.Estimate(1); beta(2) = predCurve.Coefficients.Estimate(2);
    Xo = (ES_sum - beta(1)) / beta(2);
    
    if isempty(categorical); cont = 1; else; cont = 0; end;
    
    if cont == 1 %predictor is not categorical
        [ciU, ciL, predU, predL, X, regError]  = regression_line_ci(0.05,beta,xi,yi, wi,Xo, lme, cont);
        plotRegressionConfidence(ciL, ciU, predL, predU, X, beta, xi, yi, wi,  predictors{1}, outcome{1}, properties.model, modelStats, properties, dataSet)
        if properties.ISRsep == 1
            % plot ISR results in separate figure
            
            if length(ISR)>0
                figure;
                for k = 1:length(ISR)
                    % plot ISR linear fits
                    ISRplotA(k) =  plot ( [ISR(k).([predictors{1} '_xPreda'])], [ISR(k).([predictors{1} '_yPreda'])], 'LineWidth', 1); hold on;
                    % plot observational data to which ISR was fitted
                    ISRplotB(k) =  plot ( [ISR(k).([predictors{1} '_xAdj'])], [ISR(k).([predictors{1} '_yAdj'])], 'o', 'Color',  get(ISRplotA(k), 'Color'), 'MarkerFaceColor', get(ISRplotA(k), 'Color')); hold on;
                end
            end
        else
            for k = 1:length(ISR)
                %plot ISR linear fits (overlayed on meta-regression)
                ISRplotA(k) =  plot ( [ISR(k).([predictors{1} '_xPreda'])], [ISR(k).([predictors{1} '_yPreda'])], 'r', 'LineWidth', 1); hold on;
            end
        end
    end
else
    for i = 1:length(predictors)
        exFig = figure;
        lme
        h1 = plotAdjustedResponse(lme, predictors{i}); title(['adjusted: ' predictors{i}]);
        x1=get(h1,'Xdata'); y1=get(h1,'Ydata'); close(exFig);
        
        xAdj = x1{1}; xPreda = x1{2}; yAdj = y1{1}; yPreda = y1{2};
        
        predCurve = fitlm(xPreda, yPreda);
        beta(1) = predCurve.Coefficients.Estimate(1); beta(2) = predCurve.Coefficients.Estimate(2);
        Xo = (ES_sum - beta(1)) / beta(2);
        
        if ~any(i == categorical); cont = 1; else; cont = 0; end; % check if predictor is categorical
        if cont == 1 %plot if predictor is not categorical
            
            [ciU, ciL, predU, predL, X, regError]  = regression_line_ci(0.05,beta,xAdj,yAdj, wi,Xo, lme, []);
            plotRegressionConfidence(ciL, ciU, predL, predU, X, beta, xAdj, yAdj, wi,  predictors{i}, outcome{1}, properties.model, modelStats, properties, dataSet)
            
            if properties.ISRsep == 1
                figure;
                for k = 1:length(ISR)
                    ISRplotA(k) =  plot ( [ISR(k).([predictors{i} '_xPreda'])], [ISR(k).([predictors{i} '_yPreda'])], 'LineWidth', 1); hold on;
                    ISRplotB(k) =  plot ( [ISR(k).([predictors{i} '_xAdj'])], [ISR(k).([predictors{i} '_yAdj'])], 'o', 'Color',  get(ISRplotA(k), 'Color'), 'MarkerFaceColor', get(ISRplotA(k), 'Color')); hold on;
                end
            else
                for k = 1:length(ISR)
                    ISRplotA(k) =  plot ( [ISR(k).([predictors{i} '_xPreda'])], [ISR(k).([predictors{i} '_yPreda'])], 'r',  'LineWidth', 2); hold on;
                end
                
            end
        end
    end
end

end

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
        %         es(i) = log(xr(i)) - log(xc(i));
        %         se(i) = sqrt(((sdr(i)^2)/(nr(i)*(xr(i)^2))) + ((sdc(i)^2)/(nc(i)*(xc(i)^2))));
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

function [ciU, ciL, predU, predL, X, regError] = regression_line_ci(alpha,beta,x,y,w,Xo, lme, cont)

N = length(x);

if(length(x) ~= length(y));  error(message('regression_line_ci:x and y size mismatch')); end

x_min = min(x) - abs(0.1*(max(x)-min(x)));
% if x_min > 0; x_min = 0; end;

x_max = max(x) + abs(0.1*(max(x)-min(x)));
[tW, tX] = prepareCurveData(x,w);
Xo = sum(tX .* tW)/(sum(tX)); %weighted xmean

n_pts = 200;
X = x_min:(x_max-x_min)/n_pts:x_max;
Y = ones(size(X))*beta(1) + beta(2)*X;

for i  = 1:length(x)
    [~, iMin(i)] = min(abs(x(i) - X));
end
yPred = Y(iMin);
sx = sqrt(1./w);

if length(yPred) ~= length(y) ; yPred = predict(lme); end
yResp = y;

for i = 1:length(x);
    incomplete(i) = 0;
    if isnan(yPred(i)) | isnan(yResp(i)) |...
            isempty(yPred(i)) | isempty(yResp(i)) |...
            ~isnumeric(yPred(i)) | ~isnumeric(yResp(i))
        incomplete(i) = 1;
    end
end

x = x(~incomplete); y = y(~incomplete); w = w(~incomplete); yPred = yPred(~incomplete); yResp = yResp(~incomplete);
[yResp, yPred, w, x] = prepareSurfaceData(yResp, yPred, w, x);

% SE_y_cond_x_orig = sum((y - beta(1)*ones(size(y))-beta(2)*x).^2)/(N-2); %% original error term
% SE_y_cond_x_weighted = sum(w.*((y - (beta(1)*ones(size(y))-beta(2)*x)).^2))/(sum(w)*(N-2)); %% weighted error term

SE_y_cond_x_orig = sum((yResp - yPred).^2)/(N-2); %% original error term
SE_y_cond_x_weighted = sum(w.*((yResp - yPred).^2))/(sum(w)*(N-2)); %% weighted error term
% SE_y_cond_x_monte = sum((Residual).^2)/(N-2); %% monte carlo generated error terms

SE_y_cond_x = SE_y_cond_x_orig;

%% weighted MSE
W = w/sum(w);

wMSE = sum(W.*((yResp - yPred).^2));
SSX_orig = (N-2)*var(x);
SSX = sum(((x-Xo).^2));

%% confidence calculation method
orig = 3;
if orig == 1 % original approach
    SE_Y_conf = SE_y_cond_x*(ones(size(X))*(1/N + (Xo^2)/SSX) + (X.^2 - 2*Xo*X)/SSX); % orig
    SE_Y_pred = SE_y_cond_x*(ones(size(X))*(1 + (1/N) + (Xo^2)/SSX) + (X.^2 - 2*Xo*X)/SSX); % orig
    YpredCI = (finv(1-alpha,2,N-2)*SE_Y_pred).^0.5;
    YconfCI = (finv(1-alpha,2,N-2)*SE_Y_conf).^0.5;
end;
if orig == 2; %conventioanl approach
    SE_Y_conf = SE_y_cond_x * sqrt((1/N) + (((Xo - X).^2)/((N-1)*SSX)));
    SE_Y_pred = SE_y_cond_x * sqrt(1+(1/N) + (((Xo - X).^2)/((N-1)*SSX)));
    YpredCI = tinv((1-(alpha/2)), N-2) .* SE_Y_pred;
    YconfCI = tinv((1-(alpha/2)), N-2) .* SE_Y_conf;
end
if orig == 3; % alternative approach
    SE_Y_conf = sqrt(wMSE*((1/N) + (((Xo - X).^2)./(SSX))));
    SE_Y_pred = sqrt(wMSE*(1+(1/N) + (((Xo - X).^2)./(SSX))));
    
    YpredCI = tinv((1-(alpha/2)), N-2) * SE_Y_pred;
    YconfCI = tinv((1-(alpha/2)), N-2) * SE_Y_conf;
    
    %         Ypred = (finv(1-alpha,2,N-2)*SE_Y_pred).^0.5;
    %     Yconf = (finv(1-alpha,2,N-2)*SE_Y_conf).^0.5;
    
end

ciU = Y + YconfCI;
ciL = Y - YconfCI;

predU = Y + YpredCI;
predL = Y - YpredCI;

regError.SExy_orig = SE_y_cond_x_orig;
regError.SExy_weighted = SE_y_cond_x_weighted;
regError.SSX = SSX;
regError.SEy_conf = SE_Y_conf;
regError.N = N;
end

function struc2xls(filename,s,varargin)
%STRUC2XLS Saves a data structure to an Excel file
%   Jeff Evans Dec-22-2008

%Defaults
sheet= 'Sheet1'; col= 'A'; fstrow= 1; orient= 'V';

%Optional arguments
if ~isempty(varargin)
    for j= 1:2:length(varargin)
        switch varargin{j}
            case 'Sheet'; sheet= varargin{j+1};
            case 'Col'; col= varargin{j+1};
            case 'Row'; fstrow= varargin{j+1};
            case 'Orientation'; orient= varargin{j+1};
            otherwise; error ('Unrecognized argument name');
        end
    end
end

%Output range
rangeOut=strcat(col,num2str(fstrow));
%Transform to cell
c= struct2cell(s);
%Field names
names= fieldnames(s);
%Concatenate field names and data. Create output dataset in vertical or
%horizontal orientation.
if strmatch(orient,'V')==1
    out=[names';c'];
elseif strmatch(orient,'H')==1
    out=[names, c];
end
%write
xlswrite(filename,out,sheet,rangeOut);
end