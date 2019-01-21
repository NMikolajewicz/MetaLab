function data = dataDistribution_GUI(input)

% DATADISTRIBUTION  Plot weighted histogram and normality probability plots
%                   to evaluate meta-analytic data distributions

%   [DATACOLLECTION] = DATADISTRIBUTION(DATACOLLECTION, EFFECTSIZE, DATATRANSFORMATION, TAU2ESTIMATOR, EXCLUDESTUDIES)
%   SAVESET) reads data from DATACOLLECTION and plots histograms and
%   normality probability plots of the EFFECTSIZE of interest, according
%   the DATATRANSFORMATION.
%   Optionally, alternate Tau2 estimator can be specifies.
%   Optionally, studies can be excluded from plots with EXCLUDESTUDIES.

%   Input Arguments:
%
%   DATACOLLECTION
%           String specifying name of input data colection (see PREPDATA). (ex. 'myDataCollection')
%   EFFECTSIZE
%           String specifying effect measure of interest.
%           EFFECTSIZE can be:
%           'absolute'      - evaluates distribution of absolute effect. If
%                             difference computed if negative control is
%                             reported
%           'normalized'    - evaluates distribution of normalized effect.
%                             *basal control must be reported.
%           'standardized'  - evaluates distribution of standardized effect, using Hedge's g.
%                             *basal control must be reported.
%           'ratio'         - evaluates distribution of response ratio.
%                             *basal control must be reported.
%   DATATRANSFORMATION
%           String that specifies data transformation
%           DATATRANSFORMATION can be:
%           'raw'           - no transformation
%           'log'           - log10 transformation
%   TAU2ESTIMATOR
%           String that specifies which tau2 estimator is used in for
%           random effect model
%           TAU2ESTIMATOR can be:
%           'DL'            - DerSimonian-Laird estimator (Recommended)
%           'HS'            - Hunter-Schmidt estimator
%           'H'             - Hedges estimator
%           'HM'            - Hatung-Makambi estimator
%           'SJ'            - Sidik-Jonkman estimator
%   PRECISIONMEASURE
%           String that specifies precision measure used for y-axis of
%           funnel plot
%           avaialble options:
%           'ISE'           - Inverse Standard Error (1/se)
%           'SE'            - untransformed error (plotted in reverse order)
%           'N'             - sampple size 
%   EXCLUDESTUDIES
%           Numeric array specifying which studies to remove. If all
%           studies are to be included, EXCLUDESTUDIES = [].
%   WEIGHTS
%           String specifying the precision measure used for effect size
%           calculations. options:
%           'IV'            - inverse variance (default)
%           'IVS'           - inverse standardized variance
%           'N'             - sample size
%
%   last update: 14.09.17
%---------------------------------------------------------------------------
try;
    
input

try; properties.dataCollection = input.dataCollection; catch; error('import data not specified'); end
try; properties.effectSize = input.effectSize; catch; properties.effectSize = 'absolute'; end
try; properties.dataTransformation = input.dataTransformation; catch; properties.dataTransformation = 'raw';end
try; properties.tau2estimator = input.tau2estimator; catch; properties.tau2estimator = 'DL';end
try; properties.precisionMeasure = input.precisionMeasure; catch; properties.precisionMeasure = 'ISE';end
try; properties.excludeStudies = input.excludeStudies; catch; properties.excludeStudies = [];end
try; properties.weights = input.weights; catch; properties.weights = 'IV'; end;
try; properties.nClusters = input.nClusters; catch; properties.nClusters = 2; end;
try; properties.A1 = input.A1; catch; properties.A1 = false; end
try; properties.A2 = input.A2; catch; properties.A2 = false; end
try; properties.A3 = input.A3; catch; properties.A3 = false; end
try; properties.A4 = input.A4; catch; properties.A4 = false; end
try; properties.A5 = input.A5; catch; properties.A5 = false; end
try; properties.A6 = input.A6; catch; properties.A6 = false; end
try; properties.A7 = input.A7; catch; properties.A7 = false; end
try; properties.A8 = input.A8; catch; properties.A8 = false; end

properties
load(properties.dataCollection);

for k = 1:length(D)
%     for k = 1;
    %% part 1: data extraction
    Sheet = D(k).description;                                               % dataset name
    data = dataExtraction(D(k).data, [properties.excludeStudies]);          % extact data
    covariates = D(k).covariates;                                           % extract covariates
    assignin('base', 'D', D);
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
    
    try; if properties.A1; evaluateDistribution (data, properties, Sheet); end; catch; display('error with distribution evaluation'); end                                   % evaluate weighted distributions
     try; if properties.A2; multimodalStratification(data, covariates, properties, Sheet); end; catch;display('error with covariate-cluster analysis'); end  
     try; if properties.A3; funnelPlot(data, properties, Sheet); end; catch; display('error with funnel plot');end                                                          % evaluate funnel plot
     try; if properties.A4; varianceAnalysis(data, properties, Sheet); end; catch; display('error with variance analysis');  end                                            % variance analysis
     try; if properties.A5; weightsAnalysis(data, properties, Sheet); end; catch; display('error with analysis of weights');end                                             % analysis of weighting schemes
     try; if properties.A8;  comparet2Estimator(data, covariates, properties, Sheet); end ; catch; display('error with analysis of heterogeneity estimators');end           % comparison of tau2 estimators
     try; if properties.A6; baujatAnalysis(data, properties, Sheet); end; catch; display('error with baujat plot');end                                                      % baujat plot
    %     monteCarloSampleEstimtation(data, properties, Sheet);                                                                                                             % monte carlo data resampling
     try; if properties.A7;  [exSens, cumulativeExSens] = sensitivityAnalysis(data, properties, Sheet); end; catch; display('error with sensitivity analysis');end          % sensitivty analyses
    
    
    
    
    
%       if properties.export                                                % exports results
%             fields2rmv = {'xi', 'wi_fe', 'sei_fe', 'ni','id'};              % remove array elements from data structure
%             for r = 1:length(fields2rmv);
%                 try; FINAL = rmfield(FINAL, fields2rmv{r}); catch; end;     % remove specified arrays
%                 try; exSens = rmfield(exSens, fields2rmv{r}); catch;  end
%             end
%             inputSheet = ['inputData_' SheetName{1}];                       % input sheet name
%             [inputSheet] = truncateSheet(inputSheet);                       % truncate sheetname if exceeds limits
%     
    %             eS_Sheet = ['eS_' SheetName{1}];                                
%             [eS_Sheet] = truncateSheet(eS_Sheet);
%             esCum_Sheet = ['eS_cumul_' SheetName{1}]; 
%             [esCum_Sheet] = truncateSheet(esCum_Sheet);

%             try; writetable(struct2table(exSens), properties.file, 'Sheet', eS_Sheet); catch; end
%             try; writetable(struct2table(cumulativeExSens), properties.file, 'Sheet', esCum_Sheet); catch; end
%         end



end

catch; 
      msgbox({'Error while running heterogeneity module',...
      '                                                              ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');    
end
end

function  comparet2Estimator(data, covariates, properties, Sheet);

% data transformation (if specified)
[z,sz,ID, Ni, ind] = transData (data, properties, Sheet);

% assign data
wi = 1./(sz.^2); zw = z.*wi;
S.ESi = z; S.sei_fe = sz; S.wi_fe = wi;

curInd = 1;

% store original estimator
origT2 = properties.tau2estimator;

%DL
properties.tau2estimator = 'DL';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');
curInd = curInd + 1;
 
 %HS
properties.tau2estimator = 'HS';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');
curInd = curInd + 1;

%H
properties.tau2estimator = 'H';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');
curInd = curInd + 1;

%HM
properties.t2estimator = 'HM';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');
curInd = curInd + 1;

%SJ
properties.tau2estimator = 'SJ';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');
curInd = curInd + 1;

%PM
properties.tau2estimator = 'PM';
heterogeneity = tau2Estimator(S, properties);
t2Est(curInd) = heterogeneity.t2; clear('heterogeneity');


figure; 
plot (1:length(t2Est), t2Est, 'ro','MarkerFaceColor','red');
xlabel('t^2 estimators'); ylabel('between-study variance, t^2');  title('Comparison of t^2 Estimators');
 
 set(gca, 'Xtick', 1:1:6);
set(gca,'xticklabel', {'DL', 'HS', 'H', 'HM', 'SJ', 'PM'});

 % reassign original estimator
 properties.tau2estimator = origT2;

end

function multimodalStratification(data, covariates, properties, Sheet);

% data transformation (if specified)
[z,sz,ID, Ni, ind] = transData (data, properties, Sheet);

k =properties.nClusters; % number of clusters

    
    % cluster data and assign membership
    membership = kmeans(z', k);
    
    figure; hold on;
    for i = 1:k
    histogram(z(membership == i));
    legName{i} = ['cluster '  num2str(i)];
    end
    xlabel('effect size'); ylabel('count'); legend(legName);
    title(['"' Sheet '" exploratory cluster analysis']);
    
    for i = 1:length(covariates)
        sourceID = []; curCov = []; covCodes = []; T1 = []; T2 = []; indep = [];
        unCode = []; freqC = []; rmvCode = [];
            try;
        % extract codes for given covariate
        curCov = [data.(covariates{i})];
        covCodes = curCov(ind);
        
        %pearson's chi square test for independence
        %check whether cluster membership is dependent on any of the covariates.
%         whos
        [fin(i).tbl,fin(i).chi2,fin(i).p,fin(i).labels] = crosstab(membership', covCodes);
        fin(i).covariate = covariates{i};
        
        sampleSizeCriteria(i) = sum(sum(fin(i).tbl>4));
        [n(i), m(i)] = size(fin(i).tbl);
        fin(i).sizeCriteriaSatisfied = round(100*sampleSizeCriteria(i)/(n(i)*m(i)))/100; % proportion of cells that have atleast 5 values (100% ~ valid results)
        

        unCode = unique(covCodes);
        for j = 1:length(unCode);
            freqC(j) = sum(covCodes == unCode(j));
        end
        rmvCode = unCode(freqC<=(k*4));
        
        membershipUpdated = [];
        for j = 1:length(rmvCode)
            indk = [];
            indk = find(covCodes ~= rmvCode(j));
            membershipUpdated = membership(indk);
            covCodes= covCodes(indk);
        end
        
        [fin(i).tbl2,fin(i).chi22,fin(i).p2,fin(i).labels2] = crosstab(membershipUpdated', covCodes);
        
            catch;
%                 display(['no independence test for "' covariates{i} '"']);
            end

    end
    
assignin('base', 'fin', fin);
counter1 = 1;

try;
    for i = 1:length(fin);
        if ~isnan(fin(i).p) & ~isempty(fin(i).p) & fin(i).sizeCriteriaSatisfied > 0.75
            pCov{counter1} = fin(i).covariate;
            pEl(counter1) = fin(i).p;
            counter1 = counter1 + 1;
        end
    end
    pEl = pEl(~isnan(pEl));
catch
    for i = 1:length(fin);
        if ~isnan(fin(i).p) & ~isempty(fin(i).p) 
            pCov{counter1} = fin(i).covariate;
            pEl(counter1) = fin(i).p;
            counter1 = counter1 + 1;
        end
    end
    pEl = pEl(~isnan(pEl));
end

[pEl, indCov] = sort(pEl);
pCov = pCov(indCov);
assignin('base', 'pCov', pCov);

figure;
bar(pEl, 'k'); set(gca, 'yscale', 'log'); hold on;
h1 = hline(0.05,'r'); xlabel('covariates'); ylabel('Chi-Squared Independence Test, p-value'); title(['"' Sheet '" Cluster-Covariate Dependency']);
set(gca,'xticklabel', pCov); xtickangle(45);

end

function  monteCarloSampleEstimtation(data, properties, Sheet);

% data transformation (if specified)
[z,sz,ID, Ni] = transData (data, properties, Sheet);

% heterogeneity statistisc
    S.ESi = z; S.sei_fe = sz; S.wi_fe =  1./(sz.^2); S.ni = Ni;
    heterogeneity = tau2Estimator(S, properties);
    wf = 1./(sz.^2);

% standard deviation (study-level)
sd = sz.*sqrt(Ni);

% monte-carlo sampling of study-level data
nTrials = 700;
for j = 1:nTrials
    pseudoRaw = [];
    for i = 1:length(z)
        temp = [];
        temp = z(i) + sd(i)*randn(Ni(i), 1);
        pseudoRaw = [pseudoRaw temp'];
    end
    trialMean(j) = mean(pseudoRaw);
    trialMedian(j) = median(pseudoRaw);
    trialStd(j) = std(pseudoRaw);
    overallMean(j) = mean(trialMean);
    overallStd(j) = mean(trialStd);
    overallMedian(j)  = mean(trialMedian);
    loRank(j) = (length(pseudoRaw)/2)-((1.96*sqrt(length(pseudoRaw)))/2);
    hiRank(j) = 1+(length(pseudoRaw)/2)+((1.96*sqrt(length(pseudoRaw)))/2);
    pseudoRaw = sort(pseudoRaw);
    loCi(j) = pseudoRaw(floor(loRank(j)));
    hiCi(j) = pseudoRaw(ceil(hiRank(j)));
    overallLo(j) = mean(loCi);
    overallHi(j) = mean(hiCi);
    
    
end

med.mean = overallMean(end);
med.med = overallMedian(end);
med.lo = overallLo(end);
med.hi = overallHi(end);
assignin('base', 'med', med);

%plot monte carlo simulation history
figure; 
subplot(121); plot([1:nTrials], overallMean, 'b'); xlabel('trial number'); ylabel('cumulative estimate'); hold on;  title('mean vs median convergence');
subplot(121); plot([1:nTrials], overallMedian, 'r'); xlabel('trial number');  legend('mean', 'median');
subplot(122); plot([1:nTrials], overallStd, 'k'); xlabel('trial number'); ylabel('cumulative std'); title('std convergence');

% pseudoraw and study-level effect size distributions
NM = 'probability';
figure; histogram(z, 'Normalization',NM); hold on; histogram(pseudoRaw, 'Normalization',NM); title([Sheet{1} 'Reconstructed Data']);
legend('study-level', 'reconstructed; sample-level');
xlabel('effect size'); ylabel('count');


%pseudoraw variance statistics
    % standard deviations
    pseudoStd = overallStd(end);
    pseudoSEM = pseudoStd/sqrt(length(z));
    % variance 
    pVsd = pseudoStd^2;
    pVsem = pseudoSEM^2;

%assign monte-carlo results to workspace.
monteCarlo.pseudoRaw = pseudoRaw;
monteCarlo.z = z;
monteCarlo.pseudoStd = pseudoStd;
monteCarlo.pseudoSEM = pseudoSEM;
monteCarlo.pVsd = pVsd;
monteCarlo.pVsem = pVsem;
monteCarlo.emp_t2 = heterogeneity.t2;
monteCarlo.pseudoMean = overallMean(end);

assignin('base', 'monteCarlo', monteCarlo);
end




function weightsAnalysis(data, properties, Sheet);

% log transfomration (if specified)
[z,sz,ID, Ni] = transData (data, properties, Sheet);

% assign data
wi = 1./(sz.^2); zw = z.*wi;
S.ESi = z; S.sei_fe = sz; S.wi_fe = wi;

holdThis = properties.weights;
% heterogeneity statistics (inverse variance method)
properties.weights = 'IV';
heterogeneity = tau2Estimator(S, properties);

% heterogeneity statistics (inverse standardized variance method)
properties.weights = 'IVS'
heterogeneityCV = tau2Estimator(S, properties);
properties.weights = holdThis;

% assign heterogeneity statistics
% Q = heterogeneity.Q;
t2 = heterogeneity.t2;
CV = abs(sz./z);

%fixed and random effects weights (absolute)
fixedW = 1./(sz.^2); randomW = 1./((sz.^2)+t2);

%fe, re and sample size weights (normalized)
fixedW = fixedW/(sum(fixedW));                                              % fixed effects weights
randomW = randomW/(sum(randomW));                                           % random effects weights
sizeW = Ni; sizeW = sizeW/(sum(sizeW));                                     % sample size weights
cvW = 1./(CV.^2); cvW = cvW/(sum(cvW));                                     % standardized variance weights (coefficient of variance, using standard error as precision)
unweightedW = ones(1,length(z)) .*(1/length(z));                            % equal weights (unweighted)

% estimation of fe and re model variances using standardized variance weights
F.cvE = sum(z .* cvW)/(sum(cvW));
F.cvSEi = abs(F.cvE .*CV);
cv2feW = 1./([F.cvSEi].^2);
cv2reW = 1./(([F.cvSEi].^2)+heterogeneityCV.t2);

% distribution of weights
% figure;
% subplot(321); hist(log10(fixedW)); title(['"' Sheet{1} '" fe weights']);
% subplot(322); hist(log10(randomW)); title(['"' Sheet{1} '" re weights']);
% subplot(323); hist(log10(sizeW)); title(['"' Sheet{1} '" n weights']);
% subplot(324); hist(log10(cvW)); title(['"' Sheet{1} '" cv weights']);
% subplot(325); hist(log10(cv2feW./(sum(cv2feW)))); title(['"' Sheet{1} '" cv2feW weights']);
% subplot(326); hist(log10(cv2reW./(sum(cv2reW)))); title(['"' Sheet{1} '" cv2reW weights']);

wCorel.fixedW = fixedW;
wCorel.randomW = randomW;
wCorel.sizeW = sizeW;
wCorel.fevsre = corr(fixedW, randomW, 'type', 'spearman');
wCorel.fevsn = corr(fixedW, sizeW, 'type', 'spearman');

assignin('base', 'wCorel', wCorel);

% relationship between standardized and absolute fixed and random effects weights
% figure;
% 
% subplot(121); plot ((cv2feW/sum(cv2feW)), fixedW,'ro');  hold on;
% title(['"' Sheet{1} '" cvfe vs fe']); xlabel('cv fe weights'); ylabel('fe weights');
% equality = linspace(0,1);
% plot(equality,equality,'r','LineWidth',1); 
% xlim([min((cv2feW/sum(cv2feW))) max((cv2feW/sum(cv2feW)))]); ylim([min(fixedW) max(fixedW)]);
% legend('weights', 'line of equality', 'location', 'best');
% 
% subplot(122); plot ((cv2reW/sum(cv2reW)), randomW,'bo');  hold on;
% plot(equality,equality,'r','LineWidth',1);
% title(['"' Sheet{1} '" cvre vs re']); xlabel('cv re weights'); ylabel('re weights');
% xlim([ min((cv2reW/sum(cv2reW))) max((cv2reW/sum(cv2reW)))]); ylim([min(randomW) max(randomW)]);
% legend('weights', 'line of equality', 'location', 'best');


% illustration of relative weight contribuetion
figure;
% subplot(311);
plot ([1:length(z)], unweightedW); hold on;
plot ([1:length(z)], fixedW);
plot ([1:length(z)], randomW);
plot ([1:length(z)], sizeW);
legend('unweighted', 'fe', 're', 'n'); title (['"' Sheet{1} '" comparison of relative contribution of weights']);
xlabel ('study index'); ylabel('Proportion of total weight');
% subplot(312);
% plot ([1:length(z)], unweightedW); hold on;
% plot ([1:length(z)], fixedW);
% plot ([1:length(z)], (cv2feW./(sum(cv2feW))));
% legend('unweighted', 'fe', 'cv2feW');
% subplot(313);
% plot ([1:length(z)], unweightedW); hold on;
% plot ([1:length(z)], randomW);
% plot ([1:length(z)], (cv2reW./(sum(cv2reW))));
% legend('unweighted', 're', 'cv2reW'); xlabel('study');


% effect size estimtes under various models/weightig schemes
E.fixE = sum(z .* fixedW)/(sum(fixedW));
E.randomE = sum(z .* randomW)/(sum(randomW));
E.sizeE = sum(z .* sizeW)/(sum(sizeW));
E.unweightedE = sum(z .* unweightedW)/(sum(unweightedW));
E.fixCV_E = sum(z .* cv2feW)/(sum(cv2feW));
E.randomCV_E = sum(z .* cv2reW)/(sum(cv2reW));


% monte carlo effect size
sd = sz.*sqrt(Ni);

% monte-carlo sampling of study-level data
nTrials = 700;
for j = 1:nTrials
    pseudoRaw = [];
    for i = 1:length(z)
        temp = [];
        temp = z(i) + sd(i)*randn(Ni(i), 1);
        pseudoRaw = [pseudoRaw temp'];
    end
    trialMean(j) = mean(pseudoRaw);
    trialStd(j) = std(pseudoRaw);
    trialMedian(j) = median(pseudoRaw);
    overallMean(j) = median(trialMean);
    overallStd(j) = median(trialStd);
    overallMedian(j)  = median(trialMedian);
    loRank(j) = (length(pseudoRaw)/2)-((1.96*sqrt(length(pseudoRaw)))/2);
    hiRank(j) = 1+(length(pseudoRaw)/2)+((1.96*sqrt(length(pseudoRaw)))/2);
    pseudoRaw = sort(pseudoRaw);
    loCi(j) = pseudoRaw(round(loRank(j)));
    hiCi(j) = pseudoRaw(round(hiRank(j)));
    overallLo(j) = median(loCi);
    overallHi(j) = median(hiCi);
end
    E.mc_E = overallMean(end);
    E.mc_med = overallMedian(end);
    pseudoStd = overallStd(end);
    pseudoSEM = pseudoStd/sqrt(length(z));


Wcv = 1./(abs(F.cvE .*CV).^2);

%
reW = 1./((sz.^2)+t2);
feW = 1./((sz.^2));
F.reSEi = sqrt(1./reW);
F.feSEi = sqrt(1./feW);
F.cvSEi = abs(F.cvE .*CV);

F.sizeSEi = sqrt(sum((Ni-1).*(sz.^2))/(sum(Ni)-length(Ni)));

F.unweightedSEi = std(z)/sqrt(length(Ni));

F.reSE = sqrt(1./sum(reW));
F.feSE = sqrt(1./sum(feW));
F.cv2feSE = sqrt(1./sum(cv2feW));
F.cv2reSE = sqrt(1./sum(cv2reW));
F.mcSE = pseudoSEM;
F.mcMedup = overallHi(end);
F.mcMedlo = overallLo(end);

FINAL(1).ESre = E.fixE; FINAL(1).SEre = F.feSE;
FINAL(2).ESre = E.randomE;  FINAL(2).SEre = F.reSE;
FINAL(3).ESre = E.sizeE;  FINAL(3).SEre = F.sizeSEi;
FINAL(4).ESre = E.unweightedE;  FINAL(4).SEre = F.unweightedSEi;
FINAL(5).ESre = E.fixCV_E;  FINAL(5).SEre = F.cv2feSE;
FINAL(6).ESre = E.randomCV_E;  FINAL(6).SEre = F.cv2reSE;
FINAL(7).ESre = E.mc_E;  FINAL(7).SEre = F.mcSE;
FINAL(8).ESre = E.mc_med;  FINAL(8).SEre = (overallHi(end)-overallLo(end))/2;
for i = 1:length(FINAL); FINAL(i).critValue = 1.96; FINAL(i).df = length(z)-1; end


figure; 

switch properties.dataTransformation 
    case 'log'
        
        properties.log2raw_method = 'naive';
        FINAL = log2raw(FINAL, properties);
        transformationMethod = properties.log2raw_method;
        
        assignin('base', 'FINAL', FINAL);
     errorbar([1:4],...
    [FINAL(4).ESre_raw, FINAL(1).ESre_raw, FINAL(2).ESre_raw, FINAL(3).ESre_raw],...
     [FINAL(4).ciLre_raw, FINAL(1).ciLre_raw, FINAL(2).ciLre_raw, FINAL(3).ciLre_raw],...
     [FINAL(4).ciUre_raw, FINAL(1).ciUre_raw, FINAL(2).ciUre_raw, FINAL(3).ciUre_raw],...
    'ro',...
    'MarkerFaceColor','red');   
        
        
    case'raw'
zCrit = 1.96;
errorbar([1:4],...
    [E.unweightedE, E.fixE, E.randomE, E.sizeE],...
    [zCrit*F.unweightedSEi, zCrit*F.feSE, zCrit*F.reSE, zCrit*F.sizeSEi],...
    [ zCrit*F.unweightedSEi, zCrit*F.feSE, zCrit*F.reSE, zCrit*F.sizeSEi],...
    'ro',...
    'MarkerFaceColor','red');

end
title(['"' Sheet{1} '" comparison of effect size estimators']); ylabel('effect size +- 95% CI, raw scale');
xlabel('weighting schemes');

% set(axr1, 'Ytick', -0.1:0.1:1);
% set(axr1, 'YtickLabel', {' ', '0', '0.1', '0.2', '0.3', '0.4', '0.5','0.6', '0.7', '0.8','0.9', '1'});

set(gca, 'Xtick', 1:1:4);
set(gca,'xticklabel', {'unweighted', 'FE', 'RE', 'N'});

modelEv = modelEval(E, z, F, sz, Sheet);
assignin('base', 'F', F);
assignin('base', 'E', E);
assignin('base', 'modelEv', modelEv);

end


function pooledSE = pooledVar(weights, se, es, z);
%% notes
% understand relationship between total variance = true variance + error
% variance. 

%Reliability is a function of both the total variance and the error variance. True variance is a
%population characteristic; error variance is a characteristic of the measuring instrument. 

%Comparisons of any sort can be distorted by differential reliability of variables

v1 = sum(weights);
v2 = sum(weights.^2);

% a = sum(weights.*((es-z).^2));
a = sum(weights.*(se.^2));
% b = v1/((v1^2)-v2);

% pooledSE = sqrt(a*b);
pooledSE = sqrt(a/v1);

end

function varianceAnalysis(data, properties, Sheet);

[z,sz,ID, Ni] = transData (data, properties, Sheet);
wi = 1./(sz.^2); zw = z.*wi;

CV = abs(sz./z);
assignin('base','data', data);



figure;
subplot(121); histogram(sz.*(sqrt(Ni)), 'FaceColor', 'k', 'EdgeColor', 'k'); title (['"' Sheet{1} '" Sample Variance, SD^2']); xlabel('SD^2'); ylabel('count');
subplot(122); plot (Ni, (sz.^2), 'ko'); xlabel('N'); ylabel('SE^2'); title(['"' Sheet{1} '" Effect of Sample Size on Sampling Error']);

Var = (sz.*sqrt(Ni)).^2;
a = length(ID); N = sum(Ni);
for i = 1:length(ID);
    sT(i) = ((Ni(i)-1)*(Var(i)));
    secQ(i) = (Ni(i)-1)*log(Var(i));
    Nho(i) = 1/(Ni(i)-1);
end

sTsum = sum(sT); s2p = sTsum/(N-a);
q = ((N-a)*(log(s2p))) - sum(secQ);
c = 1 + ((1/(3*(a-1)))*(sum(Nho) - (1/(N-a))));

Xo2 = (q/c);
Qinv = chi2inv(0.95, a-1);
Qp = 1-chi2cdf(Xo2, a-1);

figure;
[b, stats] = robustfit(Ni, Var);
[spearR, spearP] = corr(Ni', Var', 'type', 'Spearman');
aa.spearR = spearR; aa.spearP = spearP; aa.Ni = Ni; aa.Var = Var;
assignin('base', 'aa', aa);
plot(Ni, Var, 'ko'); hold on;
ni = linspace(min(Ni), max(Ni));
plot(ni,b(1)+b(2)*ni,'r','LineWidth',1);
Bp = [stats.p(2)];
xlabel('Sample Size (n)'); ylabel('Sample Variance'); title(['"' Sheet{1} '" Sample Size & Variance']);

if Qp < 0.05; Qtest = 1; else; Qtest = 0; end
if Bp < 0.05; Btest = 1; else; Btest = 0; end

if Qtest && Btest
    result = {'Variance dependent on sample size'; 'Study-level heteroscedasticity, Bartlet p < 0.05';'Potential problem with reporting of sample errors'; 'Data consistent with random effects';['spearman rho = ' num2str(spearR)]; ['spearman p = ' num2str(spearP)]};
elseif Qtest==1 && Btest==0
    result = {'Variance independent of sample size';'Study-level heteroscedasticity, Bartlet p < 0.05';'Data consistent  with random effects';['spearman rho = ' num2str(spearR)]; ['spearman p = ' num2str(spearP)]};
elseif Qtest==0 && Btest==1
    result = {'Variance dependent on sample size';'Study-level homoscedasticity, Bartlet p > 0.05'; 'Potential problem with reporting of sample errors';'Consider sample size-weighting scheme';['spearman rho = ' num2str(spearR)]; ['spearman p = ' num2str(spearP)]};
elseif Qtest==0 && Btest==0
    result = {'Variance independent of sample size';'Study-level homoscedasticity, Bartlet p > 0.05';'Data consistent with fixed effect';['spearman rho = ' num2str(spearR)]; ['spearman p = ' num2str(spearP)]};
end

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(2), ((ylim(2)-ylim(1))*0.75) + ylim(1),result, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');
legend('Study-level data', ['n vs var, p = ' num2str(Bp)], 'location', 'northeast');
end
%--------------------------------------------------------------------------
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
        w(i) = 1/(se(i)^2);
        esw(i) = es(i)*w(i);
        
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
            w(i) = 1/(se(i)^2);
            esw(i) = es(i)*w(i);
            
            data(i).Ni = N(i);
            data(i).SEi = se(i);
            data(i).ESi = es(i);
%             data(i).Wi = w(i);
%             data(i).ESWi = esw(i);
        catch ME
            data(i).Ni = [];
            data(i).SEi = [];
            data(i).ESi = [];
%             data(i).Wi = [];
%             data(i).ESWi = [];
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
        data(i).SDci = sdc(i);
        data(i).SDrxi = sdr(i);
        data(i).Spooli =  se(i) * sqrt(nr(i)+nc(i));
        data(i).SEi = se(i);
        data(i).ESi = es(i);
%         data(i).Wi = w(i);
%         data(i).ESWi = esw(i);
    catch ME;
    end
end
end

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
%         data(i).Wi = w(i);
    catch ME;
    end
end
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
        %         es(i) = log(xr(i)) - log(xc(i));
        %         se(i) = sqrt(((sdr(i)^2)/(nr(i)*(xr(i)^2))) + ((sdc(i)^2)/(nc(i)*(xc(i)^2))));
        es(i) = xr(i)/xc(i);
        se(i) = sqrt(((xr(i)^2)/(xc(i)^2))* (((sdr(i)^2)/(nr(i)*(xr(i)^2))) + ((sdc(i)^2)/(nc(i)*(xc(i)^2)))));
        % standard error computed using tailor approximation, under assumption cov(xr,xc) = 0.
        w(i) = 1/(se(i)^2);
        data(i).Ni = N(i);
        data(i).SEi = se(i);
        data(i).ESi = es(i);
%         data(i).Wi = w(i);
    catch;
    end
end
end



function [data] = dataExtraction(data, exStudies);
n = 1;
xr = []; e = []; ser = []; xc = []; sec = []; nc = [];
for i = 1:length(data)
    if ~isempty(exStudies); exclude = any(exStudies == data(i).ID); else; exclude = 0; end
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
            try; nc(n) = data(i).nc; catch; nc(n) = nrx(n); end
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

function funnelPlot(data, properties, Sheet);

%transform data (if log tranform selected)
[z,sz, ID, Ni] = transData (data, properties, Sheet);

% compute heterogeneity
wi = 1./(sz.^2); zw = z.*wi; S.ESi = z; S.sei_fe = sz; S.wi_fe = wi; S.ni = Ni;
heterogeneity = tau2Estimator(S, properties);

%update with new weights
sz = []; sz = [heterogeneity.SEi];
wi = []; wi = 1./(sz.^2);

switch properties.precisionMeasure
    case 'N'
        wi_re = Ni;
        S.ESre = sum(wi_re .* z)/sum(wi_re);
        S.SEre = sqrt(sum((Ni-1).*(sz.^2))/(sum(Ni)-length(Ni)));
        
        esFE = []; seFE = [];
        esRE =  S.ESre; seRE = S.SEre;
    otherwise
        %fixed effects
        ES_fix = sum(z .* wi)/(sum(wi));
        SE_fix = sqrt(1/sum(wi)); % original
        
        % random effects
        try; wi_re = 1./(heterogeneity.t2 + (sz.^2)); catch; wi_re = []; end
        sei_re = sqrt(1./wi_re);
        S.ESre = sum(wi_re .* z)/sum(wi_re);
        S.SEre = sqrt(1/sum(wi_re));
        
        % define assign to new variables
        esFE = ES_fix; seFE = SE_fix;
        esRE =  S.ESre; seRE = S.SEre;
end



% funnel plot precision measure
switch properties.precisionMeasure
    case 'SE'
        prf = sz;
        prr = sei_re;
        
        maxP = max(prf); minP = min(prf);
        % ES ref line
        Yes = linspace(minP, maxP, 100000);
        Xes = ones(1,length(Yes))*esRE;
        Xesf = ones(1,length(Yes))*esFE;
        
        % 95% CI ref lines
        XciL95 = Xes - 1.96.*(Yes); XciU95 = Xes + 1.96.*(Yes);
        
        
    case 'ISE'
        prf = 1./sz; % precision (FE);
        prr = 1./sei_re; % precision (RE);
        
        maxP = max(prf); minP = min(prf);
        % ES ref line
        Yes = linspace(minP, maxP, 100000);
        Xes = ones(1,length(Yes))*esRE;
        Xesf = ones(1,length(Yes))*esFE;
        
        % 95% CI ref lines
        XciL95 = Xesf - 1.96.*(1./Yes); XciU95 = Xesf + 1.96.*(1./Yes);
    case 'N'
        prf = Ni;
       
        maxP = max(prf); minP = min(prf);
        % ES ref line
        Yes = linspace(minP, maxP, 100000);
        Xes = ones(1,length(Yes))*esRE;
%         Xesf = ones(1,length(Yes))*esFE;
        
        % 95% CI ref lines
        XciL95 = Xes - 1.96.*(Yes); XciU95 = Xes + 1.96.*(Yes);
end

figure;

funP.z = z;
funP.se = prf;
assignin('base', 'funP', funP);
plot (z, prf,'ko'); hold on;                                                % study-level data
switch properties.precisionMeasure
    case 'N'
        plot (Xes, Yes, 'r--');  hold on;                                   % random effect (red)
    otherwise
        plot (Xesf, Yes, 'b--'); hold on;                                   % fixed effects (blue)
        plot (Xes, Yes, 'r--');  hold on;                                   % random effect (red)
        plot (XciL95, Yes, 'k--'); plot (XciU95, Yes, 'k--');               % 95% CI
end

xlim([min(z) max(z)]); ylim([0 inf]);                                       % axes limits
xlabel('Effect Size');                                                      % x axis label
title (['"' Sheet '" Funnel Plot']);                                        % title


% reverse y axis if standard error precision is used
h = gca;
switch properties.precisionMeasure
    case 'SE'
        set(h, 'yDir', 'reverse');
        ylabel('Precision (SE)');
        legend('Study-Level Data', 'Fixed Effects', 'Random Effects', 'location', 'best');              % figure legend
    case 'ISE'
        ylabel('Precision (1/SE)');
        legend('Study-Level Data','Fixed Effects', 'Random Effects', 'location', 'best');              % figure legend
    case 'N'
        ylabel('Precision (N)');
        legend('Study-Level Data','Random Effects', 'location', 'best');              % figure legend
end
end

function [z,sz,ID, Ni, ind] = transData (data, properties, Sheet)

switch properties.dataTransformation
    case 'log';
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi) &&  data(i).ESi > 0;
                sx(i) = data(i).SEi;
                ESi(i) = data(i).ESi;
                sz(i) = sqrt(log10(((sx(i)^2 )/ (ESi(i)^2)) + 1));
                z(i) = log10(ESi(i)) - (0.5 * (sz(i)^2));
                ID(i) = data(i).ID;
                Ni(i) = data(i).Ni;
                ind(i) = i;
            end
        end
    case 'raw'
        for i = 1:length(data);
            if   ~isempty(data(i).SEi) &&  ~isempty(data(i).ESi)
                sz(i) = data(i).SEi;
                z(i) = data(i).ESi;
                ID(i) = data(i).ID;
                Ni(i) = data(i).Ni;
                ind(i) = i;
            end
        end
        sz = sz(z~=0); z = z(z~=0);
end

keeper = 1;
for i = 1:length(z);
    if ~isempty(z(i)) && ~isnan(z(i)) && sz(i)~=0
        tempz(keeper) = z(i);
        tempsz(keeper) = sz(i);
        tempID(keeper) = ID(i);
        tempN(keeper) = Ni(i);
        tempInd(keeper) = ind(i);
        keeper = keeper + 1;
    end
end
z = []; sz = []; ID = []; Ni =[]; ind = [];
z = tempz; sz = tempsz; ID= tempID; Ni = tempN; ind = tempInd;
clear('tempz'); clear('tempsz'); clear('tempID'); clear('tempN'); clear('tempInd');
end

function evaluateDistribution (data, properties, Sheet);

[z,sz,ID, Ni] = transData (data, properties, Sheet);

if  ~isempty(z) & ~isempty(sz)
    figure;
    S.ESi = z; S.sei_fe = sz; S.wi_fe =  1./(sz.^2); S.ni = Ni;
    heterogeneity = tau2Estimator(S, properties);
    sz = []; sz = [heterogeneity.SEi]; 
    wf = 1./(sz.^2);
    wr = 1./((sz.^2)+heterogeneity.t2);
    
    noW = ones(1, length(wf));
    
    ewf = z.*wf;
    ewr = z.*wr;
    ewN = z.*[S.ni];
    enoW = z.*noW;

    [nBins, ~, ~] = sshist(z);    

    if ~isempty(nBins)
    [hwf, intf] = histwc(z, wf,  nBins); hwf = hwf/sum(hwf);
    [hwr, intr] = histwc(z, wr,  nBins); hwr = hwr/sum(hwr);
    [hwn, intn] = histwc(z, [S.ni],  nBins); hwn = hwn/sum(hwn);
    [hwnoW, intnoW] = histwc(z, noW,  nBins); hwnoW = hwnoW/sum(hwnoW);
    
%     subplot(3,2,[1,3,5]); 
    hold on;
    xq1 = min(intnoW):abs((max(intnoW)-min(intnoW))/1000):max(intnoW);
    s1 = interp1(intnoW,hwnoW,xq1, 'pchip');
    his1 = plot(xq1, s1, 'k');
    set(his1, 'LineWidth', 2)
    title([Sheet ' Data Distributions']); ylabel('probability density'); xlabel(['effect size, ' properties.dataTransformation]);
    xq2 = min(intf):abs((max(intf)-min(intf))/1000):max(intf);
    s2 = interp1(intf,hwf,xq2, 'pchip');
    his2 = plot(xq2, s2, 'b');
    set(his2, 'LineWidth', 2)
    xq3 = min(intr):abs((max(intr)-min(intr))/1000):max(intr);
    s3 = interp1(intr,hwr,xq3, 'pchip');
    his3 = plot(xq3, s3, 'r');
    set(his3, 'LineWidth', 2)
    
    xq4 = min(intn):abs((max(intn)-min(intn))/1000):max(intn);
    s4 = interp1(intn,hwn,xq4, 'pchip');
    his4 = plot(xq4, s4, 'g');
    set(his4, 'LineWidth', 2) 
    
    xlim([-inf inf]); ylim([0 inf]);
    legend('unweighted', 'fixed effects', 'random effects', 'sample-size');
    
%     subplot(322);      np1 = normplot(enoW);
%     xlabel('effect size'); set(np1(1), 'Marker', '.', 'MArkerEdgeColor', 'k'); title('Normal Probability Plot: Unweighted Data');
%     set(np1(2), 'LineStyle', '-', 'Color', 'k'); set(np1(3), 'LineStyle', '-', 'Color', 'k'); grid off
%     
%     subplot(323);      his2 = bar(intf, hwf, 'b');
%     ylim([0 1]);  title('Histogram: Weighted Data (Fixed)'); ylabel('weighted probability (fixed)'); xlabel('effect size'); his2.FaceAlpha = .5;
    
%     subplot(324);      np2 = normplot(ewf);
%     xlabel('effect size * weight(fixed)');set(np2(1), 'Marker', '.', 'MarkerEdgeColor', 'b'); title('Normal Probability Plot: Weighted Data (Fixed)');
%     set(np2(2), 'LineStyle', '-', 'Color', 'b'); set(np2(3), 'LineStyle', '-', 'Color', 'b'); grid off
    
%     subplot(325);      his3 = bar(intr, hwr, 'r');
%     ylim([0 1]); title('Histogram: Weighted Data (Random)'); ylabel('weighted probability (random)'); xlabel('effect size'); his3.FaceAlpha = .5;
%     subplot(326);      np3 = normplot(ewr);
%     xlabel('effect size * weight(random)'); set(np3(1), 'Marker', '.',  'MArkerEdgeColor', 'r'); title('Normal Probability Plot:  Weighted Data (Random)');
%     set(np3(2), 'LineStyle', '-', 'Color', 'r'); set(np3(3), 'LineStyle', '-', 'Color', 'r'); grid off
%     [ax,h]=subtitle(Sheet);
    
%     
%         subplot(321);      his1 = bar(intnoW, hwnoW, 'k');
%     ylim([0 1]); title('Histogram: Unweighted Data'); ylabel('probability'); xlabel('effect size'); his1.FaceAlpha = .5;
%     subplot(322);      np1 = normplot(enoW);
%     xlabel('effect size'); set(np1(1), 'Marker', '.', 'MArkerEdgeColor', 'k'); title('Normal Probability Plot: Unweighted Data');
%     set(np1(2), 'LineStyle', '-', 'Color', 'k'); set(np1(3), 'LineStyle', '-', 'Color', 'k'); grid off
%     subplot(323);      his2 = bar(intf, hwf, 'b');
%     ylim([0 1]);  title('Histogram: Weighted Data (Fixed)'); ylabel('weighted probability (fixed)'); xlabel('effect size'); his2.FaceAlpha = .5;
%     subplot(324);      np2 = normplot(ewf);
%     xlabel('effect size * weight(fixed)');set(np2(1), 'Marker', '.', 'MarkerEdgeColor', 'b'); title('Normal Probability Plot: Weighted Data (Fixed)');
%     set(np2(2), 'LineStyle', '-', 'Color', 'b'); set(np2(3), 'LineStyle', '-', 'Color', 'b'); grid off
%     subplot(325);      his3 = bar(intr, hwr, 'r');
%     ylim([0 1]); title('Histogram: Weighted Data (Random)'); ylabel('weighted probability (random)'); xlabel('effect size'); his3.FaceAlpha = .5;
%     subplot(326);      np3 = normplot(ewr);
%     xlabel('effect size * weight(random)'); set(np3(1), 'Marker', '.',  'MArkerEdgeColor', 'r'); title('Normal Probability Plot:  Weighted Data (Random)');
%     set(np3(2), 'LineStyle', '-', 'Color', 'r'); set(np3(3), 'LineStyle', '-', 'Color', 'r'); grid off
%     [ax,h]=subtitle(Sheet);
    else
        display(['note enough data in set "' Sheet{1} '" to generate distributions']);
    end
end

end

function [histw, vinterval] = histwc(vv, ww, nbins)
minV  = min(vv);
maxV  = max(vv);
delta = (maxV-minV)/nbins;
vinterval = linspace(minV, maxV, nbins)-delta/2.0;
histw = zeros(nbins, 1);
for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last' );
    if ~isempty(ind)
        histw(ind) = histw(ind) + ww(i);
    end
end
end

function [histw, vinterval] = histwcv(vv, ww, nbins)
minV  = min(vv);
maxV  = max(vv);
delta = (maxV-minV)/nbins;
vinterval = linspace(minV, maxV, nbins)-delta/2.0;
histw = zeros(nbins, 1);
indX  = arrayfun(@(xx) find(vinterval < vv(xx), 1, 'last'), 1:length(vv));
arrayfun(@(xx) evalin('caller', ['histw(indX(', sprintf('%d', xx),')) = histw(indX(', sprintf('%d', xx),')) + ww(', sprintf('%d', xx),');']), 1:length(vv));
end

function heterogeneity = tau2Estimator(S, properties)
estimator = properties.tau2estimator;

ESi = [S.ESi];
assignin('base', 'ESi', ESi);

switch properties.weights
    case 'IV'
        Wi = [S.wi_fe];
        heterogeneity.SEi = 1./sqrt([S.wi_fe]);
    case 'IVS'
        SEi = 1./sqrt([S.wi_fe]);
        assignin('base', 'SEi', SEi);
        CVi = abs(SEi./ESi);  SEi = [];
        Wi = 1./(CVi.^2);
        heterogeneity.SEi = abs(CVi * (sum((Wi.*ESi))/sum(Wi)));
        Wi = []; 
        Wi = 1./([heterogeneity.SEi].^2); 
    case 'N'
        Wi = [S.ni];
        heterogeneity.SEi = 1./sqrt([S.wi_fe]);
end


clear('S');
n = length(ESi);
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
Q = sum((Wi).*((ESi - ES_fix).^2)); %heterogeneity statistic

end

function [ax,h]=subtitle(text)
ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title(text);
if (nargout < 2)
    return
end
h=get(ax,'Title');
end

function [optN, C, N] = sshist(x,N)
% [optN, C, N] = sshist(x,N)
%
% Function `sshist' returns the optimal number of bins in a histogram
% used for density estimation.
% Optimization principle is to minimize expected L2 loss function between
% the histogram and an unknown underlying density function.
% An assumption made is merely that samples are drawn from the density
% independently each other.
%
% The optimal binwidth D* is obtained as a minimizer of the formula,
% (2K-V) / D^2,
% where K and V are mean and variance of sample counts across bins with width D.
% Optimal number of bins is given as (max(x) - min(x)) / D*.
%
% For more information, visit
% http://2000.jukuin.keio.ac.jp/shimazaki/res/histogram.html
%
% Original paper:
% Hideaki Shimazaki and Shigeru Shinomoto
% A method for selecting the bin size of a time histogram
% Neural Computation 19(6), 1503-1527, 2007
% http://dx.doi.org/10.1162/neco.2007.19.6.1503
%
% Example usage:
% optN = sshist(x); hist(x,optN);
%
% Input argument
% x:    Sample data vector.
% N (optinal):
%       A vector that specifies the number of bins to be examined.
%       The optimal number of bins is selected from the elements of N.
%       Default value is N = 2:50.
%       * Do not search binwidths smaller than a sampling resolution of data.
%
% Output argument
% optN: Optimal number of bins.
% N:    Bin numbers examined.
% C:    Cost function of N.
%
% See also SSKERNEL
%
% Copyright (c) 2009 2010, Hideaki Shimazaki All rights reserved.
% http://2000.jukuin.keio.ac.jp/shimazaki

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Setting
x = reshape(x,1,numel(x));
x_min = min(x);
x_max = max(x);

if nargin < 2
    buf = abs(diff(sort(x)));
    dx = min(buf(logical(buf ~= 0)));
    N_MIN = 2;              % Minimum number of bins (integer)
    % N_MIN must be more than 1 (N_MIN > 1).
    N_MAX = min(floor((x_max - x_min)/(2*dx)),50);
    % Maximum number of bins (integer)
    N = N_MIN:N_MAX;        % # of Bins
end

SN = 30;                    % # of partitioning positions for shift average
D = (x_max - x_min) ./ N;   % Bin Size Vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the Cost Function
Cs = zeros(length(N),SN);
for i = 1: length(N)
    
    shift = linspace(0,D(i),SN);
    for p = 1 : SN
        edges = linspace(x_min+shift(p)-D(i)/2,...
            x_max+shift(p)-D(i)/2,N(i)+1);   % Bin edges
        
        ki = histc(x,edges);               % Count # of events in bins
        ki = ki(1:end-1);
        
        k = mean(ki);                      % Mean of event count
        v = sum( (ki-k).^2 )/N(i);         % Variance of event count
        
        Cs(i,p) = ( 2*k - v ) / D(i)^2;    % The Cost Function
    end
    
end
C = mean(Cs,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Bin Size Selectioin
[Cmin idx] = min(C);
optN = N(idx);                         % Optimal number of bins
%optD = D(idx);                         % *Optimal binwidth
%edges = linspace(x_min,x_max,N(idx));  % Optimal segmentation

end



function FINAL = log2raw(FINAL, properties)
transformationMethod = properties.log2raw_method;

for i = 1:length(FINAL)
    zr = []; szr2 = []; n = [];  szr = [];
    try; zf = []; szf2 = [];  szf = []; zf = FINAL(i).ESfe; szf = FINAL(i).SEfe; szf2 = FINAL(i).SEfe^2;catch; end
    zr = FINAL(i).ESre;
    szr = FINAL(i).SEre;
    szr2 = FINAL(i).SEre^2;
    n = FINAL(i).df + 1;
    crit = FINAL(i).critValue;
    switch transformationMethod
        case 'tailor' % Higgins (2008)
            xr = 10^(zr+(szr2/2));
            varXr = (1/n)* ((10.^((2*zr)+(szr2))) * (szr2) * (1+((szr2)/2)));
            sxr = sqrt(varXr);
            xlor = xr - (crit * sxr);  xhir = xr + (crit * sxr);
            try;
                xf = 10^(zf+(szf2/2));
                varXf = (1/n)* ((10.^((2*zf)+(szf2))) * (szf2) * (1+((szf2)/2)));
                sxf = sqrt(varXf);
                xlof = xf - (crit * sxf);  xhif = xf + (crit * sxf);
            catch; end
            transMeth = 'tailor approximation';
        case 'bias' % P.Rothery (1988)
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);
            xr = (10.^zr) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            xlor = (10.^zlor) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            xhir = (10.^zhir) +  ((1.15 * szr.^2).*(n-1 ./(n)));
            sxr = (xhir-xlor)/(2*crit);
            
            try;
                zlof = zf - (crit*szf);
                zhif = zf + (crit*szf);
                xf = (10.^zf) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                xlof = (10.^zlof) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                xhif = (10.^zhif) +  ((1.15 * szf.^2).*(n-1 ./(n)));
                sxf = (xhif-xlof)/(2*crit);
            catch; disp('prob1'); end
            
            transMeth = 'bias correction';
        case 'naive'
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);
            sxr = (10^((2*zr) + szr2))* (10^(szr2) - 1);
            xr = 10^(zr+(szr2/2));
            xlor = 10^(zlor+(szr2/2)); xhir = 10^(zhir+(szr2/2));
            
            try;
                zlof = zf - (crit*szf);
                zhif = zf + (crit*szf);
                sxf = (10^((2*zf) + szf2))* (10^(szf2) - 1);
                xf = 10^(zf+(szf2/2));
                xlof = 10^(zlof+(szf2/2));xhif = 10^(zhif+(szf2/2));
            catch;  end
            transMeth = 'naive transformation';
        case 'geometric'
            zlor = zr - (crit*szr); zhir = zr + (crit*szr);
            xr = 10^(zr);
            xlor = 10^(zlor);  xhir = 10^(zhir);
            sxr = (xhir-xlor)/(2*crit);
            
            try;
                xf = 10^(zf);
                zlof = zf - (crit*szf); zhif = zf + (crit*szf);
                xlof = 10^(zlof);  xhif = 10^(zhif);
                sxf = (xhif-xlof)/(2*crit);
            catch; end
            transMeth = 'geometric mean';
    end
    
    FINAL(i).ESre_raw = xr;
    FINAL(i).SEre_raw = sxr;
    FINAL(i).ciLre_raw = xlor;
    FINAL(i).ciUre_raw = xhir;
    try;
        FINAL(i).ESfe_raw = xf;
        FINAL(i).SEfe_raw = sxf;
        FINAL(i).ciLfe_raw = xlof;
        FINAL(i).ciUfe_raw = xhif;
    catch; end
    FINAL(i).log2rawMethod = transMeth;
end

end



function baujatAnalysis(data, properties, Sheet);

% tranform data (if specified)
[z,sz,ID, Ni] = transData(data, properties, Sheet);

S.ESi = z; S.sei_fe = sz; S.wi_fe =  1./(sz.^2); S.ni = Ni;
heterogeneity = tau2Estimator(S, properties);

sz = []; sz = [heterogeneity.SEi]; 
wi = 1./(sz.^2);

ES_fix = sum(z .* wi)/(sum(wi));
Qi = ((z - ES_fix).^2) .* wi;

for i = 1:length(ID)
    ESex = []; Wiex = [];
    ESex = z; ESex(i) = [];  Wiex = wi; Wiex(i) = [];
    ESexcluded = sum(ESex .* Wiex)/(sum(Wiex));
    Zinf(i) = ((ES_fix - ESexcluded)^2)/(1/sum(Wiex));
    Qinf(i) = ((z(i)-ES_fix)^2)/(sz(i)^2);
end
%     cutoff = chi2inv(0.95, length(ID)-1);
cutoff = chi2inv(0.95, 1);
hetInd= []; hetInd = find (Qinf>cutoff);
hetStudies = ID(hetInd);
percentHet = 100*length(hetStudies)/length(Qi);

figure;
for i = 1:length(Qinf)
    if Qinf(i) >= cutoff
        plot(Qinf(i), Zinf(i), 'ro'); hold on;
    else
        plot(Qinf(i), Zinf(i), 'ko'); hold on;
    end
end

% for i = 1:length(hetInd)
%     text(Qinf(hetInd(i)), Zinf(hetInd(i)),num2str(ID(hetInd(i))),  'HorizontalAlignment', 'left');
% end

% switch properties.dataTransformation
%     case 'log'
%         set(gca, 'yscale', 'log');
% end

% Th = round(100*(sum(Qinf >= cutoff) / length(Qinf)));
% h = vline(cutoff,'r',[num2str(Th) '% studies above 95th percentile {\chi}^2']);

xlabel('Contribution to Overall Heterogeneity'); ylabel('Influence on Overal Result'); title (['"' Sheet{1} '" Baujat Plot']);
end


function hhh=vline(x,in1,in2)
%
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
    end
    
    
    
    
    g=ishold(gca);
    hold on
    
    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end
    
    if g==0
        hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else

if nargout
    hhh=h;
end
end




function [totalExSens, cumulativeExSense] = sensitivityAnalysis(data, properties, Sheet);

% transform data (if specified)
[z,sz, ID, Ni] = transData (data, properties, Sheet);

% heterogeneity statistics

S.ESi = z; S.sei_fe = sz; S.wi_fe =  1./(sz.^2); S.ni = Ni;
heterogeneity = tau2Estimator(S, properties);

%  fixed effects estimate
sz = []; sz = abs([heterogeneity.SEi]);
wi = 1./(sz.^2); zw = z.*wi;
ES_fix = sum(z .* wi)/(sum(wi));
SE_fix = sqrt(1/sum(wi)); % original


% reassign data to "total" structure
for i = 1:length(ID)
    total(i).id = ID(i);
    total(i).sei_fe = sz(i);
    total(i).wi_fe = wi(i);
    total(i).xi = z(i);
    total(i).ni = Ni(i);
end
idNum = [total.id];

SEi = [total.sei_fe];
Wi =  [total.wi_fe];
ESi = [total.xi];  Ni = [total.ni];

for k = 1:length(ESi)-2
    exSens = [];
    for i = 1:length(idNum)
        ESiEx = []; WiEx = [];   NiEx = []; seEx = []; sdEx = [];
        ESiEx =ESi; WiEx = Wi ;   NiEx = Ni; seEx = 1./(sqrt(WiEx)); sdEx = seEx;
        ESiEx(i)=[]; WiEx(i)=[] ; NiEx(i)=[]; seEx(i)=[]; sdEx(i)=[];
        subgroup = [];
        subgroup.y = ESiEx; subgroup.se = seEx; subgroup.n = NiEx; subgroup.indicatorNumber = idNum(i);
        S.subgroup = idNum(i);
        S= dataSummary(subgroup, S);
        S = EffectSizeEstimator(S, properties);
        S.excludeStudyID = S.subgroup;
        S = rmfield(S,'subgroup');
        try; exSens = [exSens S]; catch; exSens = S; end
    end
    Qi = []; Qmin = []; Imin = [];
    
    Qi = [exSens.Q];
    [~, Imin] = min(Qi);
    
    
    
    hetStudy(k) = exSens(Imin).excludeStudyID;
    hetStudyQ(k) = exSens(Imin).Q;
    hetStudyI(k) = exSens(Imin).I2;
    hetStudyES(k) = exSens(Imin).ESre;
    hetStudySE(k) = exSens(Imin).SEre;
    hetStudycritValue(k) = exSens(Imin).critValue;
    hetStudyhetP(k) = exSens(Imin).Qp;
    hetStudyt2(k) = exSens(Imin).t2;
    hetStudyciWidth(k) = exSens(Imin).ciWidth_re;
    hetStudydf(k) = length(ESiEx)-1;
    hetStudyH2(k) = hetStudyQ(k)/ hetStudydf(k);
    
    ESi(Imin) = []; Wi(Imin) = []; Ni(Imin)=[]; idNum(Imin) = [];
    
    if k == 1; totalExSens = exSens;
        fields2rmv = {'ESi', 'wi_fe', 'ni',  'id'};
        totalExSens = rmfield(exSens, fields2rmv);
        
        for j = 1:length(totalExSens)
            totalExSens(j).H2 = totalExSens(j).Q / (length(ESi)-1);
        end
    end
    
end

for i = 1:length(hetStudy)
    cumulativeExSense(i).ExcludedStudyID = hetStudy(i);
    cumulativeExSense(i).Q = hetStudyQ(i);
    cumulativeExSense(i).I2 = hetStudyI(i);
    cumulativeExSense(i).ESre = hetStudyES(i);
    cumulativeExSense(i).SEre = hetStudySE(i);
    cumulativeExSense(i).critValue = hetStudycritValue(i);
    cumulativeExSense(i).ciWidth_Re = hetStudyciWidth(i);
    cumulativeExSense(i).Qp = hetStudyhetP(i);
    cumulativeExSense(i).t2 = hetStudyt2(i);
    cumulativeExSense(i).df = hetStudydf(i);
    cumulativeExSense(i).H2 = hetStudyH2(i);
end

%% prep data for graphical representation
ssy = [totalExSens.ESre];
sss = [totalExSens.SEre];
cr = [totalExSens.critValue];
ciW = sss.*cr;
ssx = [1:length(ssy)];
ssp = [totalExSens.Qp];

ES_cum = [ cumulativeExSense.ESre];
ciWidth = [ cumulativeExSense.SEre] .* [cumulativeExSense.critValue ];
cumStudy = [1:length(ES_cum)];
Qp = [cumulativeExSense.Qp];

Th = round(100*(sum(Qp < 0.05) / length(Qp)));

%% heterogeneity sensitivity
figure;
subplot(121);
yyaxis left; plot(ssx, [totalExSens.H2], 'r');  ylabel('H^2', 'Color', 'r'); axl1 = gca; axl1.YColor = 'r';
yyaxis right; plot(ssx, [totalExSens.I2], 'k'); ylabel('I^2', 'Color', 'k'); axl1 = gca; axl1.YColor = 'k';
xlim([min(ssx), max(ssx)]);
xlabel('single studies excluded');
title([Sheet 'Single Study Exclusion Analysis']);
subplot(122);
yyaxis left; plot(cumStudy, [cumulativeExSense.H2], 'r');  ylabel('H^2', 'Color', 'r'); axl1 = gca; axl1.YColor = 'r';
yyaxis right; plot(cumStudy, [cumulativeExSense.I2], 'k'); ylabel('I^2', 'Color', 'k'); axl1 = gca; axl1.YColor = 'k';
xlim([-inf inf]);
xlabel('cumulative studies excluded');
title([Sheet ' Cumulative Exclusion Analysis']);

%% effect size sensitivity
figure;
subplot(121);
yyaxis left
[~,~]=jbfill_Reg(ssx,(ssy-ciW),(ssy+ciW),[1 0 0],[1 0 0],0,0.3); hold on;
ylabel('effect size', 'Color', 'r');
axl1 = gca; axl1.YColor = 'r';

yyaxis right
p1 = plot(ssx, ssp, 'k', 'LineWidth', 1);
axr1 = gca; axr1.YColor = 'k';
ylabel('Q test, p-value', 'Color', 'k');
% set(axr1, 'yscale', 'log');
% assignin('base', 'totalExSens', totalExSens);
ylim([-0.1 1]);
set(axr1, 'Ytick', -0.1:0.1:1);
set(axr1, 'YtickLabel', {' ', '0', '0.1', '0.2', '0.3', '0.4', '0.5','0.6', '0.7', '0.8','0.9', '1'});
% xlim([-inf inf]);
xlim([min(ssx), max(ssx)]);
xlabel('single studies excluded');
title([Sheet 'Single Study Exclusion Analysis']);

subplot(122);


yyaxis left
[~,~]=jbfill_Reg(cumStudy,(ES_cum-ciWidth),(ES_cum+ciWidth),[1 0 0],[1 0 0],0,0.3); hold on;
ylabel('effect size', 'Color', 'r');
axl2 = gca; axl2.YColor = 'r';
yyaxis right
p2 = plot(cumStudy, Qp, 'k', 'LineWidth', 1);
axr2 = gca; axr2.YColor = 'k';
% set(axr2, 'yscale', 'log');
ylabel('Q test, p-value', 'Color', 'k');


h = vline(sum(Qp < 0.05),'k',['T_h =' num2str(Th) '%']);

ylim([0 1]);
% xlim([-inf inf]);
xlim([-inf, inf]);
xlabel('cumulative studies excluded');
title([Sheet 'Cumulative Exclusion Analysis']);


end


function S = dataSummary(subgroup, S)
ES = [subgroup.y]; SEM = [subgroup.se]; Ni = [subgroup.n]; indicatorNum = [subgroup.indicatorNumber];

df = length(ES) - 1;
if df > 0
    W_fix = 1./(SEM.^2);
    ESW = ES.*W_fix;
    
    ES_fix = sum(ESW) / sum(W_fix);
    Q = sum((W_fix).*((ES - ES_fix).^2));
    
    H2 = Q / df; %higgins & thompson 2002, AKA Birge's ratio (Birge 1932).
    if Q > df;  I = 100*(H2 - 1)/H2;
    else; I = 0; end
    
    if Q > df+1; SE_lnH = (0.5)*((log(Q)-log(df))/((sqrt(2*Q))-(sqrt((2*(df+1))-3))));
    else; SE_lnH = sqrt((1/(2*(df-1)))*(1-(1/(3*((df-1)^2))))); end
    
    H2_ciU = exp((log(H2) + (1.96 * SE_lnH)));
    H2_ciL = exp((log(H2) - (1.96 * SE_lnH)));
    I_ciU =  real(100*(H2_ciU - 1)/H2_ciU);
    I_ciL =  real(100*(H2_ciL - 1)/H2_ciL);
    if I_ciU < 0; I_ciU = 0; end; if isnan(I_ciU); I_ciU = 0; end
    if I_ciL < 0; I_ciL = 0; end; if isnan(I_ciL); I_ciL = 0; end
    
else;
    try; W_fix = 1./(SEM.^2);I = []; I_ciL = []; I_ciU = []; Q = []; H2 = [];
    catch; I = []; I_ciL = []; I_ciU = []; Q = []; W_fix = []; H2 = [];  end
end


% p_Q = 1-(round(1000*chi2cdf(Q,df))/1000);
p_Q = 1-(chi2cdf(Q,df));

S.I2 = I;
S.I2lo = I_ciL;
S.I2hi = I_ciU;
S.H2 = H2;
S.df = df;
S.Q = Q;
S.Qp = p_Q;
S.ESi = ES;
S.wi_fe = W_fix;
S.sei_fe = SEM;
S.ni = Ni;
S.id = indicatorNum;

end


function [S] = EffectSizeEstimator(S, properties)

if S.df ~=0
    try heterogeneity = tau2Estimator(S, properties); catch ME
        getReport(ME)
    end
    
    xi = [S.ESi];
    sei_fe = [S.sei_fe];
    wi_fe = [S.wi_fe];
    k = length(xi);
    df = k-1;
    t2 = heterogeneity.t2;
    
    ES_fix = sum(xi .* wi_fe)/(sum(wi_fe));
    SE_fix = sqrt(1/sum(wi_fe));
    S.C = sum(wi_fe) - ((sum((wi_fe).^2))/sum(wi_fe));
    
    try; wi_re = 1./(heterogeneity.t2 + (sei_fe.^2)); catch; wi_re = []; end
    sei_re = sqrt(1./wi_re);
    
    S.ESre = sum(wi_re .* xi)/sum(wi_re);
    S.SEre = sqrt(1/sum(wi_re));
    S.ESfe = ES_fix;
    S.SEfe = SE_fix;
    
    properties.CIestimator = 'z'; %update this to be included in user prompt window (291017)
    switch properties.CIestimator
        
        case 'z'
            z = norminv(1-0.05/2) ;
            S.ciWidth_re = z * sqrt(1/sum(wi_re));
            S.ciWidth_fe =  z * sqrt(1/sum(wi_fe));
            S.critValue = z;
            S.ciEstimator = 'zdist';
            
        case 't'
            t = tinv(1-0.05/2, df);
            S.ciWidth_re = t * sqrt(1/sum(wi_re));
            S.ciWidth_fe = t * sqrt(1/sum(wi_fe));
            S.critValue = t;
            S.ciEstimator = 'tdist';
            
        case 'QA'
            b = 2.061 + (4.902/k) + (0.756/sqrt(k)) - (0.959/log(k));
            S.ciWidth_re = b * sqrt(1/sum(wi_re));
            S.ciWidth_fe = b * sqrt(1/sum(wi_fe));
            S.critValue = b;
            S.ciEstimator = 'quantileApproximation';
    end
    
    S.tauEstimator = properties.tau2estimator;
    S.t2 = t2;
    
else
    S.C = [];
    S.ESre = [S.xi];
    S.SEre = [S.sei_fe];
    S.ESfe = [S.xi];
    S.SEfe = [S.sei_fe];
    switch properties.CIestimator
        case 'z'
            z = norminv(1-0.05/2) ;
            S.ciWidth_re = z *  S.SEre;
            S.ciWidth_fe =  z *  S.SEre;
            S.critValue = z;
            S.ciEstimator = 'zdist';
            
        case 't'
            t = tinv(1-0.05/2, df);
            S.ciWidth_re = t *  S.SEre;
            S.ciWidth_fe = t *  S.SEre;
            S.critValue = t;
            S.ciEstimator = 'tdist';
            
        case 'QA'
            b = 2.061 + (4.902/k) + (0.756/sqrt(k)) - (0.959/log(k));
            S.ciWidth_re = b *  S.SEre;
            S.ciWidth_fe = b *  S.SEre;
            S.critValue = b;
            S.ciEstimator = 'quantileApproximation';
    end
    S.tauEstimator = 'not applicable';
    S.t2 = [];
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

function modelEv = modelEval(E, z, F, sz, Sheet);

fields = fieldnames(E);

%evaluate bias
for i = 1:length(fields)
   modelEv.([fields{i} '_bias']) = mean(E.(fields{i}) - z);
   a(i) =  modelEv.([fields{i} '_bias']);
end

 
% evaluate MSE
for i = 1:length(fields)
   modelEv.([fields{i} '_MSE']) = mean((E.(fields{i}) - z).^2);
   b(i) =  modelEv.([fields{i} '_MSE']);
end

% evaluate coverage

ciLi = z-(1.96*sz); ciUi = z+(1.96*sz);

    % assign data

    lo_feE = E.fixE - 1.96*F.feSE; hi_feE = E.fixE + 1.96*F.feSE; 
    lo_reE = E.randomE - 1.96*F.reSE; hi_reE = E.randomE + 1.96*F.reSE;
    lo_feCvE = E.fixCV_E - 1.96*F.cv2feSE; hi_feCvE = E.fixCV_E + 1.96*F.cv2feSE;
    lo_reCvE = E.randomCV_E - 1.96*F.cv2reSE; hi_reCvE = E.randomCV_E + 1.96*F.cv2reSE;
    lo_nE = E.sizeE - 1.96*F.sizeSEi; hi_nE = E.sizeE + 1.96*F.sizeSEi;
    lo_mcE = E.mc_E - 1.96*F.mcSE; hi_mcE = E.mc_E + 1.96*F.mcSE;
    lo_mcMed = F.mcMedlo; hi_mcMed = F.mcMedup;
    lo_noWE = E.unweightedE - 1.96*F.unweightedSEi;   hi_noWE = E.unweightedE + 1.96*F.unweightedSEi;

    
    % compute coverage according to interval overlap
    covType = 1;
    if covType == 1; % coverage means confidence intervals overlap
        for j = 1:length (z)
            cover_fe(j) = 1;
            cover_re(j) = 1;
            cover_feCV(j) = 1;
            cover_reCV(j) = 1;
            cover_n(j) = 1;
            cover_mc(j) = 1;
            cover_mcMed(j) = 1;
            cover_noweight(j) = 1;
            
            if lo_feE > ciUi(j);
                cover_fe(j) = 0; end;
            if ciLi(j) > hi_feE;
                cover_fe(j) = 0; end;
            
            if lo_reE > ciUi(j); cover_re(j) = 0; end;
            if ciLi(j) > hi_reE; cover_re(j) = 0; end;
            
            if lo_feCvE > ciUi(j); cover_feCV(j) = 0; end;
            if ciLi(j) > hi_feCvE; cover_feCV(j) = 0; end;
            
            if lo_reCvE > ciUi(j); cover_reCV(j) = 0; end;
            if ciLi(j) > hi_reCvE; cover_reCV(j) = 0; end;
            
            if lo_nE > ciUi(j); cover_n(j) = 0; end;
            if ciLi(j) > hi_nE; cover_n(j) = 0; end;
            
            if lo_mcE > ciUi(j); cover_mc(j) = 0; end;
            if ciLi(j) > hi_mcE; cover_mc(j) = 0; end;
            
            if lo_mcMed > ciUi(j); cover_mcMed(j) = 0; end;
            if ciLi(j) > hi_mcMed; cover_mcMed(j) = 0; end;
            
            if lo_noWE > ciUi(j); cover_noweight(j) = 0; end;
            if ciLi(j) > hi_noWE; cover_noweight(j) = 0; end;
        end
    elseif covType ==2 ; % coverage means estimate confidence intervals cover study level means
        for j = 1:length (z)
            cover_fe(j) = 1;
            cover_re(j) = 1;
            cover_feCV(j) = 1;
            cover_reCV(j) = 1;
            cover_n(j) = 1;
            cover_mc(j) = 1;
            cover_mcMed(j) = 1;
            cover_noweight(j) = 1;
            
            if lo_feE > z(j);
                cover_fe(j) = 0; end;
            if z(j) > hi_feE;
                cover_fe(j) = 0; end;
            
            if lo_reE > z(j); cover_re(j) = 0; end;
            if z(j) > hi_reE; cover_re(j) = 0; end;
            
            if lo_feCvE > z(j); cover_feCV(j) = 0; end;
            if z(j) > hi_feCvE; cover_feCV(j) = 0; end;
            
            if lo_reCvE > z(j); cover_reCV(j) = 0; end;
            if z(j) > hi_reCvE; cover_reCV(j) = 0; end;
            
            if lo_nE > z(j); cover_n(j) = 0; end;
            if z(j) > hi_nE; cover_n(j) = 0; end;
            
            if lo_mcE > z(j); cover_mc(j) = 0; end;
            if z(j) > hi_mcE; cover_mc(j) = 0; end;
            
            if lo_mcMed > z(j); cover_mcMed(j) = 0; end;
            if z(j) > hi_mcMed; cover_mcMed(j) = 0; end;
            
            if lo_noWE > z(j); cover_noweight(j) = 0; end;
            if z(j) > hi_noWE; cover_noweight(j) = 0; end;
        end
    end
    
    % store coverage results
    coverage(1) = 100*sum(cover_noweight)/length(z);
    coverage(2) = 100*sum(cover_fe)/length(z);
    coverage(3) = 100*sum(cover_re)/length(z);
%     coverage(3) = 100*sum(cover_feCV)/length(z);
%     coverage(4) = 100*sum(cover_reCV)/length(z);
    coverage(4) = 100*sum(cover_n)/length(z);
%     coverage(6) = 100*sum(cover_mc)/length(z);
%     coverage(7) = 100*sum(cover_mcMed)/length(z);
   


figure; 
bar(coverage); title(['"' Sheet{1} '" coverage']); 
set(gca,'xticklabel', {'unweight','FE', 'RE', 'n'});
ylim([0 100]);
% figure; 
% subplot(121); bar(a); title(['"' Sheet{1} '" bias']); 
% set(gca,'xticklabel', fields);
% subplot(122); bar(b); title(['"' Sheet{1} '" MSE']);
% set(gca,'xticklabel', fields);

end

function hhh=hline(y,in1,in2)
% function h=hline(y, linetype, label)
% 
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    x=get(gca,'xlim');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
        else
            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
        end
    end

    if g==0
    hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end
end