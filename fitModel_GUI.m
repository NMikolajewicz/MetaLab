function fitModel_GUI(input)

% FITMODEL  Fit aggregate data with model of interest (linear, quadratic, hyperbolic or hill).

%   [FITRESULTS] = FITMODEL(FILE, MODEL, EFFECTSIZE, NTRIALS, HISTORY, XBOUNDS, YBOUNDS, SAVE)
%   fits MODEL to data imported from FILE using Monte-Carlo error propagation
%   method and estimates model parameters in terms of EFFECTSIZE. Monte-Carlo
%   data are sampled NTRIALS times and restricted by limits imposed by
%   XBOUNDS and YBOUNDS.
%   Optionally, simulation HISTORY is plotted to evaluate convergence of
%   parameter estimates.
%
%   Input Arguments:
%-----------------------
%   FILE
%           String that specifies name of the Microsoft Excel Spreadsheet
%           from which data is imported. Import format is same as exported from
%           MULTIPLE export format from DATAEXTRACTION (see DATAEXTRACTION).
%           Data in spreadsheet must have following headers with each row
%           representing an independent observation:
%           'exposure'      - x exposure data
%           'xr'            - y response data
%           'ser'           - 'xr' standard error (interested in sampling error,
%                             not natural variation of 'xr', hence standard deviation
%                             is not appropriate here)
%           'nr'            - sample size of 'xr'
%
%           Optionally, basal control data can also be included. Absolute
%           EFFECTSIZE can be computed without basal data, however for
%           relative and hedges EFFECTSIZE, basal data must be provided as
%           input:
%           'xc'            - y control (basal/baseline) data
%           'sec'           - 'xc' standard error
%           'nc'            - sample size of 'xc'
%   MODEL
%           String specifying model to fit data. MODEL can be:
%           'linear'        - linear model: (b(1).*x)+b(2)
%                             minimum input observations to fit model: 3
%           'quadratic'     - quadratic model (vertex form): ((b(1).*((x-b(2)).^2)) + b(3))
%                             minimum input observations to fit model: 4
%           'expo'          - exponential model: b(1)*exp(b(2)*x)
%                             minimum input observations to fit model: 3
%           'hyperbolic     - hyperbolic function: (b(1) * x) ./ (b(2) + x)
%                             minimum input observations to fit model: 3
%           'hill'          - hill function: (b(1) * ((x.^b(3)))./ ((b(2).^b(3)) + (x.^b(3))))
%                             minimum input observations to fit model: 4
%   EFFECTSIZE
%           String specifying effect measure of interest.
%           EFFECTSIZE can be:
%           'absolute'      - evaluates model parameters on absolute scale.
%                             Absolute difference computed if negative control is reported
%           'normalized     - evaluates model parameters on normalized scale.
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%           'standardized'   - evaluates model parameters on standardized scale.
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%           'ratio'         - evaluates model parameters on ratiometric scale
%                             *basal control ('xc' 'sec', 'nc') must be
%                             included in FILE as headers
%   NTRIALS
%           Integer specifying number of Monte-Carlo samplings made to
%           estimate model parameter distriution (>250 recommended).
%           As value of NTRIALS tends to infinity, parameter estimates will
%           converge to expected values. However, as NTRIALS increases,
%           computational time increases.
%   HISTORY
%           Logical (TRUE/FALSE). Specifies whether parameter estimate
%           history is plotted. Produces plot illustrating parameter
%           estimates as function of simulation number.
%   XBOUNDS
%           Specifies range of X values included in the monte-carlo sampled
%           estimates. ex.[0 inf] will exclude negatively sampled X values
%           from being included in the Monte-Carlo estimate calculation.
%   YBOUNDS
%           Specifies range of Y values included in the monte-carlo sampled
%           estimates. ex.[-inf inf] will include all values sampled in the
%           Monte-Carlo estimate calculation.
%   EXPORT
%           Logical (TRUE/FALSE). Specifies whether fit results are
%           exported to FILE
%   ESTIMATECI
%           Logical (TRUE/FALSE). Specifies whether confidence intervals
%           are estimate via bootstrap method for fitted regression curvge
%   CIMETHOD
%           String specifying type of confidence interval estimattion. Can
%           be:
%           'sem'           - bootstraps standard error of fitted
%                             monte-carlo simulated models and calculates confidence 
%                             intervals according to z-distribution. 
%           'percentile'    - same idea as 'sem', except 2.5% and 97.5%
%                             percentiles are bootstrapped.
%
%   updated 15.09.17
% input

try;
warning('off','all');

try; properties.file = input.file; catch; error('input file not specified'); end
try; properties.model = input.model; catch;properties.model = 'linear'; end
try; properties.effectSize = input.effectSize; catch; properties.effectSize = 'absolute'; end
try; properties.nTrials = input.nTrials; catch; properties.nTrials = 500; end
try; properties.history = input.history; catch; properties.history = false; end
try; properties.yBounds = input.yBounds; catch; properties.yBounds = [-inf inf]; end
try; properties.export = input.export; catch; properties.export = false; end
try; properties.ciMethod = input.ciMethod; catch; properties.ciMethod = 'sem'; end
try; properties.estimateCI = input.estimateCI; catch; properties.estimateCI = false; end


properties

properties.xBounds = [-inf inf];
%--------------------------------------------------------------------------
properties.reqHeaders =  {'exposure', 'xr', 'ser', 'nr'};
properties.addHeaders = {'xc', 'sec', 'nc'};

[checklist, eligData] = fitEligibility(properties);

for i = 1:length(eligData)
    
    restructData = restructureData(eligData(i));
    StudyName = eligData(i).Study; sheetName = eligData(i).sheetName;
    
    switch properties.effectSize
        case 'absolute'
            data =  AbsoluteDifference(restructData);                               % absolute difference
        case 'normalized'
            data =  NormalizedDifference(restructData);                             % normalized difference
        case 'standardized'
            data =  hedgesG(restructData);
        case 'ratio'
            data = respRatio(restructData);
    end
    
    %% prep fitting variables
    [x,y,n,sdy,sdx,varY,varX, wX, wY] = prepInput(data, sheetName);
    
    %% determine initial parameter estimates
    
    [x0, lb, ub] = initialEstimations(x,y, properties.model);
    
    %% monte carlo error propagation
    [fitResults] = monteCarloErrorPropagation (x, y, varX, varY, wX, wY, x0, lb, ub,  n, properties, sheetName);
    
    estimateQuality = menu({'Assess Quality of Parameters', 'Save figures (if needed for later reference)'}, 'estimates acceptable', 'estimates not adequate');
    parameterEstimates(i) = prep2save(fitResults, StudyName, sheetName, properties, estimateQuality);
    close all; clear fitResult
end

properties.file = [properties.file '.xlsx'];
if properties.export
    SheetName = 'FitResults';
    try; writetable(struct2table(parameterEstimates), properties.file, 'Sheet',SheetName);
    catch; writetable(struct2table(parameterEstimates, 'AsArray', true), properties.file, 'Sheet', SheetName); end
    display(['Parameter estimates succesfully exported to "' properties.file '"']);
else
    disp('Fitting complete. Results were not exported');
end

clearvars('-except', 'parameterEstimates', 'checklist', 'properties','input');


catch;
    
 
  msgbox({'Error while running fit data module',...
      '                                                              ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0, lb, ub] = initialEstimations(x,y,model)

[ymax, ind] = max(y);
midX = (x(round(ind/2)));

% define model to fit
switch model
    case 'hill'
        modelFun =@(b,x)(b(1) * ((x.^b(3)))./ ((b(2).^b(3)) + (x.^b(3))));
        tryInit = [ymax midX 1];
    case 'hyper'
        modelFun = @(b,x)(((b(1) * x) ./ (b(2) + x)) + b(3));
        tryInit = [ymax midX];
    case 'expo'
        modelFun = @(b,x)(b(1).*exp((b(2)).*x));
        tryInit = [ymax midX];
    case 'quadratic'
        %         modelFun = @(b,x)(b(1).*(x.^2))+(b(2).*x) + (b(3));
        modelFun = @(b,x) ((b(1).*((x-b(2)).^2)) + b(3));
        tryInit = [0 midX 0];
    case 'linear'
        modelFun = @(b,x)(b(1).*x)+b(2);
        tryInit = [0 0];
    case 'biphasic'
        modelFun =@(b,x) b(1) + ((b(2) - b(1))./(1+ (10.^((x-b(3)).*b(4))))) + ((b(5)-b(1)) ./ (1+(10.^((b(6)-x).*b(7)))));
        %b1 = amin, b2 = amax1, b3 = x1, b4 = h1, b5 = amax2, b6 = x2, b7 = h2
        tryInit = [min(y) median(y) 0 1 max(y) 0 1];
end

% attempt initial model fit
try
    try
        [beta,~,~,~,~,~] = nlinfit(x, y,modelFun, tryInit );
    catch ME;
        getReport(ME)
        [beta,~,~,~,~,~] = nlinfit(x(2:end), y(2:end),modelFun, tryInit);
    end
catch ME; getReport(ME); display('failed to estimate initial set of parameters'); beta = []; end
alpha = beta;

% assign initial parameter estimates and parameter constraints
switch model
    case 'hill'
        if isempty(alpha)
            alpha = [ymax midX 1];
            lb = [-inf -inf -inf];
            ub = [inf inf inf];
        else;
            lb = [-inf -inf -inf];
            ub = [inf inf inf];
        end
        x0 = [alpha(1) alpha(2) alpha(3)];
    case 'hyper'
        if isempty(alpha);
            alpha = [ymax midX 0];
            lb = [-inf -inf (min(x)-(abs(5*min(x))))];
            ub = [inf inf (min(x)+(abs(5*min(x))))];
        else;
            lb = [-inf -inf (min(x)-(abs(5*min(x))))];
            ub = [inf inf (min(x)+(abs(5*min(x))))];
        end
        x0 = [alpha(1) alpha(2) alpha(3)];
    case 'expo'
        if isempty(alpha);
            alpha = [0 0];
            lb = [-inf -inf];
            ub = [inf inf ];
        else;
            lb = [-inf -inf ];
            ub = [inf inf ];
        end
        x0 = [alpha(1) alpha(2)];
    case 'quadratic'
        if isempty(alpha);
            alpha = [0 0 0];
            lb = [-inf -inf -inf];
            ub = [inf inf inf];
        else;
            lb = [-inf -inf -inf];
            ub = [inf inf inf];
        end
        x0 = [alpha(1) alpha(2)  alpha(3)];
    case 'linear'
        if isempty(alpha);
            alpha = [0 0];
            lb = [-inf -inf];
            ub = [inf inf ];
        else;
            lb = [-inf -inf ];
            ub = [inf inf ];
        end
        x0 = [alpha(1) alpha(2)];
    case 'biphasic'
        if isempty(alpha);
            alpha = [min(y) median(y) median(x) 0.5 median(y) median(x) 0.5];
            lb = [-inf -inf -inf -inf -inf -inf -inf];
            ub = [inf inf inf inf inf inf inf];
%             lb = [-100 -5 1 0 0 20 0.01];
%             ub = [50 0 10 0.5 100 60 0.06];
            % temp
            %             lb = [-100 -5 2 0 0 20 0.01];
            %             ub = [50 0 10 0.5 100 60 0.05];
            % ok
            %                         lb = [-100 -5 0 0 0 20 0.005];
            %             ub = [50 0 10 0.1 100 60 0.1];
        else;
              lb = [-inf -inf -inf -inf -inf -inf -inf];
            ub = [inf inf inf inf inf inf inf];
%             lb = [-100 -5 2 0 0 20 0.01];
%             ub = [50 0 10 0.2 100 60 0.05];
        end
        if alpha(1) > alpha(2); alpha(1) = min(y); alpha(2) = median(y);  end
        if alpha(1)  > alpha(5); alpha(1) = min(y); alpha(5) = max(y); end
        if alpha(3) < min(x); alpha(3) = median(x); end;
        if alpha(6) < min(x); alpha(6) = median(x); end;
        
        x0 = [alpha(1) alpha(2) alpha(3) alpha(4) alpha(5) alpha(6) alpha(7)];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitResults] = monteCarloErrorPropagation (x, y, vX, vY, wX, wY, x0, lb, ub, n, properties, Sheet);


nTrials = properties.nTrials;
xBound = [properties.xBounds];
yBound = [properties.yBounds];
history = properties.history;
model = properties.model;

if isempty(n); nObs = 1; else; nObs = n; end            %if sample size is undefined, assign value of 1
if length(n) == 1; nObs = ones(size(y)) * n; end;

allY = []; allX = [];
c2 = 1; c1 = 1; c3 = 1; c4 = 1; c5 = 1; c6 = 1; c7 = 1;
xq = linspace(min(x), max(x), 200);
MCrealizations = figure; subplot(121);
% for i = 1:nTrials
i = 1;
while i < nTrials+1
    try;
    setY = []; setX = []; yGen = []; xGen = [];
    for j = 1:length(y)
        try;
        tempY = y(j) + sqrt(vY(j))*randn(nObs(j), 1);
        catch; tempY = y(j); end
        yGen(j) = mean(tempY);
    end
    for j = 1:length(x)
        try;
            tempX = x(j) + sqrt(vX(j))*randn(nObs(j), 1);
        catch; tempX = x(j); end
        xGen(j) = mean(tempX);
    end
    
    l1 = (yGen > yBound(1)); l2 = (yGen < yBound(2));
    l3 = (xGen > xBound(1)); l4 = (xGen < xBound(2));
    
    for k = 1:length(yGen);
        if l1(k) == 0 | l2(k) ==0 | l3(k)==0 | l4(k)==0
            remo(k) = 1;
        else; remo(k) = 0; end
    end
    yGen = yGen(~remo); xGen = xGen(~remo);                 %remove values that exceed  x and y boundaries
    
    percntiles = prctile(yGen,[5 95]);
    outlierIndex = yGen < percntiles(1) | yGen > percntiles(2);
    yGen(outlierIndex) = [];  xGen(outlierIndex) = [];
    
    setY = yGen'; setX = xGen';
    allY = [allY; setY]; allX = [allX; setX];
    
    switch model
        case 'hill'
            modelFun = @(b)((b(1) * ((setX.^b(3)))./ ((b(2).^b(3)) + (setX.^b(3))))-setY);
            backX = [0 0 0];
        case 'hyper'
            modelFun = @(b)((((b(1) * setX) ./ (b(2) + setX))+(b(3)))-setY);
            backX = [0 0 0];
        case 'expo'
            modelFun = @(b)((b(1).*exp((b(2).*setX)))-setY);
            backX = [0 0];
        case 'quadratic'
            modelFun = @(b)(((b(1).*((setX-b(2)).^2)) + b(3))-setY); % vertex form
            backX = [0 0 0];
%             backX = [200 50 150];
        case 'linear'
            modelFun = @(b)(((b(1).*setX)+b(2))-setY);
            backX = [0 0];
        case 'biphasic'
            modelFun = @(b) ((b(1) + ((b(2) - b(1))./(1+ (10.^((setX-b(3)).*b(4))))) + ((b(5)-b(1)) ./ (1+(10.^((b(6)-setX).*b(7))))))-setY);
            backX = [0 0 0 0 0 0 0];
    end
    
    switch model
        case 'hill'
            modelFun2 =@(b,x)(b(1) * ((x.^b(3)))./ ((b(2).^b(3)) + (x.^b(3))));
        case 'hyper'
            modelFun2 = @(b,x)(((b(1) * x) ./ (b(2) + x))+b(3));
        case 'expo'
            modelFun2 = @(b,x)(b(1).*exp((b(2)).*x));
        case 'quadratic'
            modelFun2 = @(b,x) ((b(1).*((x-b(2)).^2)) + b(3)); %vertex form
        case 'linear'
            modelFun2 = @(b,x)(b(1).*x)+b(2);
        case 'biphasic'
            modelFun2 =@(b,x)(b(1) + ((b(2) - b(1))./(1+ (10.^((x-b(3)).*b(4))))) + ((b(5)-b(1)) ./ (1+(10.^((b(6)-x).*b(7))))));
    end
    
    %fit to monte carlo sampled data
    
    assignin('base', 'x0', x0);
    assignin('base', 'lb', lb);
    assignin('base', 'ub', ub);
    assignin('base', 'modelFun', modelFun);
    
    properties.fitAlgorithm = 'fitnlm';
    switch properties.fitAlgorithm
        case 'lsq'
            options = optimoptions('lsqnonlin','Display','off');
            R= [];
            try;     [beta{i},~,R,~,~,~,J]  = lsqnonlin(modelFun, x0,  lb, ub, options);
            catch;  [beta{i},~,R,~,~,~,J]  = lsqnonlin(modelFun, backX,  lb, ub, options); end
        case 'nlinfit'
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'bisquare';
            CovB = []; MSE = []; R= []; J = [];
            try;     [beta{i},R{i},J{i}, CovB, MSE]  = nlinfit(setX, setY, modelFun2, x0,  opts);
            catch;  [beta{i},R{i},J{i}, CovB, MSE]  = nlinfit(setX, setY, modelFun2, backX, opts); end
        case 'fitnlm'
%             opts = statset('nlinfit');
%             opts.RobustWgtFun = 'bisquare';
%             
            for retry = 1:50 %sometime models return nan/inf. retry fit few times to omit these errors
                try;
%                     fitR = fitnlm(setX, setY, modelFun2, x0, 'Options', opts);
                    fitR = fitnlm(setX, setY, modelFun2, x0);
                catch;
                    try;
                        assignin('base', 'modelFun2', modelFun2);
                        assignin('base', 'backX', backX);
%                         fitR = fitnlm(setX, setY, modelFun2, rand().*x0, 'Options',opts);
                         fitR = fitnlm(setX, setY, modelFun2, rand().*x0);
                    catch;
                        try;
                        backX = rand().*x0;
%                         fitR = fitnlm(setX, setY, modelFun2, backX, 'Options',opts);
                        fitR = fitnlm(setX, setY, modelFun2, backX);
                        catch; end
                    end
                end;
                try; fitR; 
                    assignin('base', 'fitR', fitR);
                     beta{i} = fitR.Coefficients.Estimate(:);
                     [~,predCI{i}] = predict(fitR, xq', 'Simultaneous', true);
                    break; catch; end;  %break out of loop if  solution is found
            end

            

           
            
            
    end
    assignin('base', 'beta', beta);
%     assignin('base', 'predCI', predCI);
    Xpred{i} = setX;
    
    %current progress
    if rem(i*10, nTrials) == 0;
        display (['current progress: ' num2str(round(1000*i / nTrials)/10)]);
    end
    
    if isreal(beta{i}(2)); beta2(c2) = beta{i}(2); c2 = c2 + 1;  else; end;
    if isreal(beta{i}(1)); beta1(c1) = beta{i}(1); c1 = c1 + 1;  else; end;
    try; if isreal(beta{i}(3)); beta3(c3) = beta{i}(3); c3 = c3 + 1;  else; end;  catch ME; getReport(ME);
    end
    try; if isreal(beta{i}(4)); beta4(c4) = beta{i}(4); c4 = c4 + 1;  else; end;  catch ME; getReport(ME);
    end
    try; if isreal(beta{i}(5)); beta5(c5) = beta{i}(5); c5 = c5 + 1;  else; end;  catch ME; getReport(ME);
    end
    try; if isreal(beta{i}(6)); beta6(c6) = beta{i}(6); c6 = c6 + 1;  else; end;  catch ME; getReport(ME);
    end
    try; if isreal(beta{i}(7)); beta7(c7) = beta{i}(7); c7 = c7 + 1;  else; end;  catch ME; getReport(ME);
    end
    
    % fit model and plot curreniteration
    modelResp = [];
    switch model
        case 'hill'
            plot (setX, setY,'o'); hold on;
            modelResp = (( beta{i}(1) .* (([1:max(x)].^beta{i}(3))))./ ((beta{i}(2).^beta{i}(3)) + ([1:max(x)].^beta{i}(3))));
            modelObs = (( beta{i}(1) .* ((xq.^beta{i}(3))))./ ((beta{i}(2).^beta{i}(3)) + (xq.^beta{i}(3))));
           
            plot([1:max(x)],modelResp);
        case 'hyper'
            plot (setX, setY,'o'); hold on;
            modelResp = (( beta{i}(1) .* ([1:max(x)]))./ (beta{i}(2) + ([1:max(x)]))) + (beta{i}(3));
            modelObs = (( beta{i}(1) .* (xq))./ (beta{i}(2) + (xq))) + (beta{i}(3));
            plot([1:max(x)],modelResp);
        case 'expo'
            plot(setX,setY,'o'); hold on;
            modelResp = (beta{i}(1).*exp((beta{i}(2)).*[1:max(x)]));
            modelObs = (beta{i}(1).*exp((beta{i}(2)).*xq));
            plot([1:max(x)],modelResp);
        case 'quadratic'
            plot (setX, setY,'o'); hold on;
            %             plot([1:max(x)],(beta{i}(1).*([1:max(x)].^2))+(beta{i}(2).*[1:max(x)]) + beta{i}(3)); %standard form
            modelResp = ((beta{i}(1).*(([1:max(x)]-beta{i}(2)).^2)) + beta{i}(3));
            modelObs = ((beta{i}(1).*((xq-beta{i}(2)).^2)) + beta{i}(3));
            plot([1:max(x)],modelResp); %vertex form
        case 'linear'
            plot (setX, setY,'o'); hold on;
            modelResp = (beta{i}(1).*[1:max(x)])+beta{i}(2);
            modelObs = (beta{i}(1).*xq)+beta{i}(2);
            plot([1:max(x)],modelResp);
        case 'biphasic'
            plot (setX, setY,'o'); hold on;
            modelResp = beta{i}(1) + ((beta{i}(2) - beta{i}(1))./(1+ (10.^(([1:max(x)]-beta{i}(3)).*beta{i}(4))))) + ((beta{i}(5)-beta{i}(1)) ./ (1+(10.^((beta{i}(6)-[1:max(x)]).*beta{i}(7)))));
            modelObs = beta{i}(1) + ((beta{i}(2) - beta{i}(1))./(1+ (10.^((xq-beta{i}(3)).*beta{i}(4))))) + ((beta{i}(5)-beta{i}(1)) ./ (1+(10.^((beta{i}(6)-xq).*beta{i}(7)))));
            plot([1:max(x)], modelResp );
    end
    xlabel('x'); ylabel('y'); title('Monte Carlo Simulations (total)');
    
    modR{i} = modelObs; %store for bootstrapping of model intervals
    
    
    % if final iteration, find and remove outlying parameter estimates
    if i == nTrials
        [~, keepInd] = rmvOutlier (beta1);
        [~, keepInd] = rmvOutlier (beta2, keepInd);
        try;[~, keepInd] = rmvOutlier (beta3,keepInd);catch; end
        try;[~, keepInd] = rmvOutlier (beta4,keepInd);catch; end
        try;[~, keepInd] = rmvOutlier (beta5,keepInd);catch; end
        try;[~, keepInd] = rmvOutlier (beta6,keepInd);catch; end
        try;[~, keepInd] = rmvOutlier (beta7,keepInd);catch; end
        assignin('base', 'keepInd', keepInd);
        beta1 = beta1(keepInd==1);
        beta2 = beta2(keepInd==1);
        try;beta3 = beta3(keepInd==1);catch; end
        try;beta4 = beta4(keepInd==1);catch; end
        try;beta5 = beta5(keepInd==1);catch; end
        try;beta6 = beta6(keepInd==1);catch; end
        try;beta7 = beta7(keepInd==1);catch; end
        
    end
        
% store mean parameter estimates
            beta1Mean(c1) = mean(beta1); beta1Std(c1) = std(beta1);
            beta2Mean(c2) = mean(beta2); beta2Std(c2) = std(beta2);
            
            try; beta3Mean(c3) = mean(beta3); beta3Std(c3) = std(beta3);   catch ME; getReport(ME);
            end
            try; beta4Mean(c4) = mean(beta4); beta4Std(c4) = std(beta4);   catch ME; getReport(ME);
            end
            try; beta5Mean(c5) = mean(beta5); beta5Std(c5) = std(beta5);  catch ME; getReport(ME);
            end
            try; beta6Mean(c6) = mean(beta6); beta6Std(c6) = std(beta6);   catch ME; getReport(ME);
            end
            try; beta7Mean(c7) = mean(beta7); beta7Std(c7) = std(beta7);  catch ME; getReport(ME);
            end

        try; clear('fitR'); end
        i = i+1;
    catch; 
        try; clear('fitR'); end
    end
end
% fitting ends here

% summary of results starts here
figure('NumberTitle', 'off', 'Name', 'Parameter Distributions');
switch model
    case 'hill'
        subplot(131); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(132); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(133); h3 = histogram(beta3); title('beta 3'); xlabel('beta 3'); ylabel('count'); set(h3, 'FaceColor', 'b', 'EdgeColor', 'b');
        x0 = [beta1, beta2, beta3];
    case 'hyper'
        subplot(131); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'k', 'EdgeColor', 'k');
        subplot(132); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'k', 'EdgeColor', 'k');
        subplot(133); h3 = histogram(beta3); title('beta 3'); xlabel('beta 3'); ylabel('count'); set(h3, 'FaceColor', 'k', 'EdgeColor', 'k');
        x0 = [beta1, beta2, beta3];
    case 'expo'
        subplot(121); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'k', 'EdgeColor', 'k');
        subplot(122); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'k', 'EdgeColor', 'k');
        x0 = [beta1, beta2];
    case 'quadratic'
        subplot(131); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'r', 'EdgeColor', 'r');
        subplot(132); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'r', 'EdgeColor', 'r');
        subplot(133); h3 = histogram(beta3); title('beta 3'); xlabel('beta 3'); ylabel('count'); set(h3, 'FaceColor', 'r', 'EdgeColor', 'r');
        x0 = [beta1, beta2, beta3];
    case 'linear'
        subplot(121); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'g', 'EdgeColor', 'g');
        subplot(122); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'g', 'EdgeColor', 'g');
        x0 = [beta1, beta2];
    case 'biphasic'
        subplot(331); h1 = histogram(beta1); title('beta 1'); xlabel('beta 1'); ylabel('count'); set(h1, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(332); h2 = histogram(beta2); title('beta 2'); xlabel('beta 2'); ylabel('count'); set(h2, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(333); h3 = histogram(beta3); title('beta 3'); xlabel('beta 3'); ylabel('count'); set(h3, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(334); h4 = histogram(beta4); title('beta 4'); xlabel('beta 4'); ylabel('count'); set(h4, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(335); h5 = histogram(beta5); title('beta 5'); xlabel('beta 5'); ylabel('count'); set(h5, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(336); h6 = histogram(beta6); title('beta 6'); xlabel('beta 6'); ylabel('count'); set(h6, 'FaceColor', 'b', 'EdgeColor', 'b');
        subplot(337); h7 = histogram(beta7); title('beta 7'); xlabel('beta 7'); ylabel('count'); set(h7, 'FaceColor', 'b', 'EdgeColor', 'b');
        x0 = [beta1, beta2, beta3, beta4, beta5, beta6, beta7];
end


%% plot history;
% if history
%     switch model
%         case 'hill'
%             figure; subplot(131); plot ([1:c1], beta1Mean,'b'); ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(132); plot ([1:c2], beta2Mean,'b'); ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end))]);
%             subplot(133); plot ([1:c3], beta3Mean,'b'); ylabel ('beta3 mean'); xlabel('Simulation Number');
%             title (['beta 3: ' num2str(beta3Mean(end))]);
%         case 'hyper'
%             figure; subplot(131); plot ([1:c1], beta1Mean,'k');ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(132); plot ([1:c2], beta2Mean,'k');ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end))]);
%             subplot(133); plot ([1:c3], beta3Mean,'k');ylabel ('beta3 mean'); xlabel('Simulation Number');
%             title (['beta 3: ' num2str(beta3Mean(end))]);
%         case 'expo'
%             figure; subplot(121); plot ([1:c1], beta1Mean,'m');ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(122); plot ([1:c2], beta2Mean,'m');ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end))]);
%         case 'quadratic'
%             figure; subplot(131); plot ([1:c1], beta1Mean,'r');ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(132); plot ([1:c2], beta2Mean,'r');ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end)) ]);
%             subplot(133); plot ([1:c3], beta3Mean,'r');ylabel ('beta3 mean'); xlabel('Simulation Number');
%             title (['beta 3: ' num2str(beta3Mean(end))]);
%         case 'linear'
%             figure; subplot(121); plot ([1:c1], beta1Mean,'g'); ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(122); plot ([1:c2], beta2Mean,'g'); ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end))]);
%         case 'biphasic'
%             figure;
%             subplot(331); plot ([1:c1], beta1Mean,'b'); ylabel ('beta1 mean'); xlabel('Simulation Number');
%             title (['beta 1: ' num2str(beta1Mean(end))]);
%             subplot(332); plot ([1:c2], beta2Mean,'b'); ylabel ('beta2 mean'); xlabel('Simulation Number');
%             title (['beta 2: ' num2str(beta2Mean(end))]);
%             subplot(333); plot ([1:c3], beta3Mean,'b'); ylabel ('beta3 mean'); xlabel('Simulation Number');
%             title (['beta 3: ' num2str(beta3Mean(end))]);
%             subplot(334); plot ([1:c4], beta4Mean,'b'); ylabel ('beta4 mean'); xlabel('Simulation Number');
%             title (['beta 4: ' num2str(beta4Mean(end))]);
%             subplot(335); plot ([1:c5], beta5Mean,'b'); ylabel ('beta5 mean'); xlabel('Simulation Number');
%             title (['beta 5: ' num2str(beta5Mean(end))]);
%             subplot(336); plot ([1:c6], beta6Mean,'b'); ylabel ('beta6 mean'); xlabel('Simulation Number');
%             title (['beta 6: ' num2str(beta6Mean(end))]);
%             subplot(337); plot ([1:c7], beta7Mean,'b'); ylabel ('beta7 mean'); xlabel('Simulation Number');
%             title (['beta 7: ' num2str(beta7Mean(end))]);
%     end
% end

betaMean(1) = beta1Mean(end); betaMean(2) = beta2Mean(end);
betaStd(1) = beta1Std(end); betaStd(2) = beta2Std(end);
try; betaMean(3) = beta3Mean(end); betaStd(3) = beta3Std(end);catch; end
try; betaMean(4) = beta4Mean(end); betaStd(4) = beta4Std(end);catch; end
try; betaMean(5) = beta5Mean(end); betaStd(5) = beta5Std(end);catch; end
try; betaMean(6) = beta6Mean(end); betaStd(6) = beta6Std(end);catch; end
try; betaMean(7) = beta7Mean(end); betaStd(7) = beta7Std(end);catch; end

keepThis = 1;
for i = 1:length(keepInd)
    if keepInd(i) == 1
    modelpredCI_LO(:,keepThis) = predCI{i}(:,1);
    modelpredCI_HI(:,keepThis) = predCI{i}(:,2);
    modelResp_q(:,keepThis) = modR{i}(:);
    keepThis = keepThis+1;
    end
end

assignin('base', 'modelResp_q', modelResp_q);
assignin('base', 'modR', modR);
assignin('base', 'setX', setX);
assignin('base', 'x', x);

figure(MCrealizations); subplot(122); plot(xq, modelResp_q);
xlabel('x'); ylabel('y'); title('Monte Carlo Simulations (outliers removed)');
 
if nTrials < 1000; bootN = 1000; else; bootN = 2*nTrials; end;

switch properties.ciMethod
    case 'percentiles'
        pctl_lo = bootstrp(bootN, @(u) [prctile(u,2.5)], modelResp_q');
        pctl_hi = bootstrp(bootN, @(u) [prctile(u,97.5)], modelResp_q');
        assignin('base', 'pctl_lo', pctl_lo);
        assignin('base', 'pctl_hi', pctl_hi);
        ciLO_Q = median(pctl_lo, 1); ciLO_Q = smooth(ciLO_Q);
        ciHI_Q = median(pctl_hi, 1);ciHI_Q = smooth(ciHI_Q);
    case 'sem'
        bootstat = bootstrp(bootN, @std, modelResp_q');
        bootsq = bootstat.^2;
        bootsum = sum(bootsq,1)/size(bootsq,1);
        bootSE = sqrt(bootsum);
        bootCI = bootSE*1.96;
        assignin('base', 'bootstat', bootstat);
        assignin('base', 'bootCI', bootCI);
end
% 




figure;
switch model
    case 'hill'
        finalY = (( betaMean(1) .* ((xq.^betaMean(3))))./ ((betaMean(2).^betaMean(3)) + (xq.^betaMean(3))));
        modelY = (( betaMean(1) .* ((x.^betaMean(3))))./ ((betaMean(2).^betaMean(3)) + (x.^betaMean(3))));
        
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        if properties.estimateCI; ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3); end
        plot(xq,finalY,'r');
        modelFunction = 'b_1(x^{b3}) / (b_2^{b3} + x^{b3})'; modelName = 'sigmoidal';
    case 'hyper'
        finalY = ((betaMean(1) .* (xq))./ (betaMean(2) + (xq))) + betaMean(3);
        modelY = ((betaMean(1) .* (x))./ (betaMean(2) + (x))) + betaMean(3);
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        if properties.estimateCI; ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3); end
        plot(xq,finalY, 'r');
        modelFunction = '(b_1(x) / (b_2 + x))+b_3'; modelName = 'hyperbolic';
    case 'expo'
        finalY = (betaMean(1).*exp((betaMean(2)).*xq));
        modelY = (betaMean(1).*exp((betaMean(2)).*x));
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3)
        plot (xq,finalY, 'r');
        modelFunction = 'b_1exp(b_2x)'; modelName = 'exponential';
    case 'quadratic'
        finalY = ((betaMean(1).*((xq-betaMean(2)).^2)) +betaMean(3));
        modelY = ((betaMean(1).*((x-betaMean(2)).^2)) +betaMean(3));
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        if properties.estimateCI; ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3); end
        plot(xq,finalY, 'r'); %vertex form
        modelFunction = 'b_1(x-b_2)^2 + b_3'; modelName = 'quadratic';
    case 'linear'
        finalY = (betaMean(1).*xq)+betaMean(2);
        modelY = (betaMean(1).*x)+betaMean(2);
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        if properties.estimateCI; ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3); end
        plot(xq,finalY, 'r');
        modelFunction = 'b_1(x) + b_2 '; modelName = 'linear';
    case 'biphasic'
        finalY = (betaMean(1) + ((betaMean(2) - betaMean(1))./(1+ (10.^((xq-betaMean(3)).*betaMean(4))))) + ((betaMean(5)-betaMean(1)) ./ (1+(10.^((betaMean(6)-xq).*betaMean(7))))));
        modelY = (betaMean(1) + ((betaMean(2) - betaMean(1))./(1+ (10.^((x-betaMean(3)).*betaMean(4))))) + ((betaMean(5)-betaMean(1)) ./ (1+(10.^((betaMean(6)-x).*betaMean(7))))));
        try; ciLO_Q = finalY -  bootCI; ciHI_Q = finalY +  bootCI; catch; end
        if properties.estimateCI; ciplot(ciLO_Q,ciHI_Q,xq,'r'); hold on;  alpha(0.3); end
        plot(xq,finalY ,'r');
        modelFunction = 'b_1(x^{b3}) / (b_2^{b3} + x^{b3}) + b_4(x^{b6}) / (b_5^{b6} + x^{b6})'; modelName = 'biphasic';
end

wplot = markerWeighting(wY);
hold on;

for i = 1:length(x)
    if isnan(wplot(i)) | isinf(wplot(i))
        wplot(i) = mean(wplot, 'omitnan');
    end
    plot(x(i), y(i), 'o', 'MarkerSize', wplot(i), 'MarkerEdgeColor', 'k');
end;
xlabel('x'); ylabel('y'); title(['final model fit: ' Sheet]);

[r2, ~] = rsquare(y,modelY);

%store fit results
fitResults.modelName = modelName;
fitResults.modelFunction = modelFunction;
fitResults.beta = betaMean;
fitResults.betaSEM = betaStd;
fitResults.r2 = r2;

str = {'FIT RESULTS', ['modelName: ' modelName], ['modelFunction: ' modelFunction], ['beta means: ' num2str(betaMean)], ['beta se: ' num2str(betaStd)], ['r^2 = ' num2str(r2)]};
TextLocation(str,'location', 'best')

end


function [x,y,n,sdy,sdx,varY,varX, wX, wY] = prepInput(data, Sheet);

%% prep fitting variables
x = [data.exposure];
y = [data.ESi];
n = [data.nr];
sdy = [data.SEi] .* n;
varY = [data.SEi] .^2;
wY = 1./([data.SEi].^2);
try;
    sdx = [data.sx] .* n;
    varX = [data.sx] .^2;
    wX = 1./([FINAL.sx].^2);
catch;
    %     display('Variability in exposure not included');
    sdx =0; varX = 0; wX  = ones(size(wY));
end

figure; errorbar (x,y, [data.SEi], 'ko');
title(['Input Data: ' Sheet]); xlabel('x'); ylabel('y');

if isrow(y); y = y'; end; if isrow(x); x = x'; end;
end

function wplot = markerWeighting(wi)
weighting = 1; %1: weighting illustrated according to FE, 2: weighting illustrated according to RE
if weighting == 1; markerScale = 50; else; markerScale = 20; end
wi(isinf(wi)) = max(wi(~isinf(wi)));
wplot = wi + (0.1* max(wi));
wplot = markerScale .* wplot./max(wplot);
end

function minObserv = modelPredictors(model)

switch model
    case 'hill'
        minObserv = 4;
    case 'hyper'
        minObserv = 4;
    case 'expo'
        minObserv = 3;
    case 'quadratic'
        minObserv = 4;
    case 'linear'
        minObserv = 3;
    case 'biphasic'
        minObserv = 8;
end

end

function [checklist, eligData] = fitEligibility(properties)

[~,availableSheets] = xlsfinfo(properties.file);
reqHeaders = properties.reqHeaders;
addHeaders = properties.addHeaders;
counter1 = 1; counter2 = 1;
for h = 1:length(availableSheets);
    [~,~, raw] = xlsread(properties.file, availableSheets{h});
    inputHeaders{h} = raw(1,:);
    checklist(h).Study = [];
    checklist(h).sheetName = availableSheets{h};
    checklist(h).outcomePresent = 0;
    checklist(h).basalPresent = 0;
    checklist(h).absElig = 0;
    checklist(h).relElig = 0;
    checklist(h).exposure = [];
    checklist(h).xr = [];
    checklist(h).ser = [];
    checklist(h).nr = [];
    checklist(h).xc = [];
    checklist(h).sec = [];
    checklist(h).nc = [];
    
    n = 1;
    % extract required headings
    for i = 1:length(reqHeaders)
        for j = 1:length(inputHeaders{h})
            if strcmp(reqHeaders{i}, inputHeaders{h}{j});
                checklist(h).outcomePresent = checklist(h).outcomePresent + 1;
                checklist(h).(reqHeaders{i}) = [raw{2:end,j}];
            end
        end
    end
    
    %extract study name if available
    for j = 1:length(inputHeaders{h})
        if strcmp('Study', inputHeaders{h}{j});
            checklist(h).Study = [raw{2,j}];
        end
    end
    
    %extract control data if present
    for i = 1:length(addHeaders)
        for j = 1:length(inputHeaders{h})
            if strcmp(addHeaders{i}, inputHeaders{h}{j});
                %              {'xc', 'sec', 'nc'};
                assignin('base', 'checklist', checklist);
                checklist(h).(addHeaders{i}) = [raw{2:end,j}];
                %                 checklist(h).(addHeaders{i})
                %                 ~isempty(checklist(h).(addHeaders{i}))
                %                 ~isnan(checklist(h).(addHeaders{i}))
                if ~isempty(checklist(h).(addHeaders{i})) & ~isnan(checklist(h).(addHeaders{i}));
                    checklist(h).basalPresent = checklist(h).basalPresent + 1;
                else
                    checklist(h).(addHeaders{i}) = [];
                end
            end
        end
    end
    
    minObserv = modelPredictors(properties.model);
    checklist(h).nPredictors = length(checklist(h).xr);
    if checklist(h).nPredictors >= minObserv & checklist(h).outcomePresent > 3
        checklist(h).absElig = 1;
        absData(counter1) = checklist(h);
        counter1 = counter1 + 1;
    end
    if checklist(h).nPredictors >= minObserv & checklist(h).basalPresent > 2
        checklist(h).relElig = 1;
        relData(counter2) = checklist(h);
        counter2 = counter2 + 1;
    end
end

absEligibility = sum([checklist.absElig]);
relEligibility = sum([checklist.relElig]);

fields2rmv = {'outcomePresent','basalPresent', 'eligible', 'nPredictors'};
for r = 1:length(fields2rmv);
    try; absData = rmfield(absData, fields2rmv{r}); catch;  end;
    try; absData = rmfield(relData, fields2rmv{r}); catch;  end;
end

switch properties.effectSize
    case 'absolute'
        eligData = absData;
    case 'normalized'
        eligData = relData;
    case 'standardized'
        eligData = relData;
    case 'ratio'
        eligData = relData;
end
end


%% absolute difference
function data =  AbsoluteDifference(data)

for i = 1:length(data);
    if ~isempty(data(i).nc) &&   ~isempty(data(i).sec) &&  ~isempty(data(i).xc)
        
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        
        if isempty(nc(i)); nc(i) = nan(); end
        
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
            es(i) = xr(i); % absolute difference
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
        
        if isempty(nc(i)); nc(i) = nan(); end
        
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
        
        if isempty(nc(i)); nc(i) = nan(); end
        
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


function output = restructureData(input);

exposure = [input.exposure];
xr = [input.xr];
ser = [input.ser];
nr = [input.nr];
xc  = [input.xc];
sec = [input.sec];
nc = [input.nc];

for i = 1:length(xr);
    try; output(i).xr = xr(i); catch; output(i).xr = []; end;
    try; output(i).ser = ser(i); catch; output(i).ser = [];  end;
    try; output(i).nr = nr(i); catch; output(i).nr = []; end;
    try; output(i).xc = xc(i); catch; output(i).xc = [];  end;
    try; output(i).sec = sec(i); catch; output(i).sec = [];  end;
    try; output(i).nc = nc(i); catch;  output(i).nc = []; end;
    try; output(i).exposure = exposure(i); catch;  output(i).exposure = []; end;
    
    if isnan(output(i).xc) | isempty(output(i).xc); output(i).xc = []; end
    if isnan(output(i).sec) | isempty(output(i).sec); output(i).sec = []; end
    if isnan(output(i).nc) | isempty(output(i).nc); output(i).nc = []; end
end


end


function hOut = TextLocation(textString,varargin)

l = legend(textString,varargin{:});
t = annotation('textbox');
t.String = textString;
t.Position = l.Position;
delete(l);
t.LineStyle = 'None';

if nargout
    hOut = t;
end
end

function parameterEstimates = prep2save(fitResults, StudyName, sheetName, properties, estimateQuality)

if estimateQuality == 1; good = 1; else; good = 0; end

parameterEstimates.Study = StudyName;
parameterEstimates.dataSheet = sheetName;
parameterEstimates.effectSize = properties.effectSize;
parameterEstimates.nSimulations = properties.nTrials;
parameterEstimates.modelName = fitResults.modelName;
parameterEstimates.modelFunction = fitResults.modelFunction;
parameterEstimates.r2 = fitResults.r2;
beta = [fitResults.beta];
betaSEM = [fitResults.betaSEM];
for i = 1:length(beta);
    parameterEstimates.(['beta' num2str(i)]) = beta(i);
end
for i = 1:length(beta);
    parameterEstimates.(['betaSE' num2str(i)]) = betaSEM(i);
end
parameterEstimates.acceptableEstimates = good;
end


%% response ratio
function data = respRatio(data);
for i = 1:length(data);
    try;
        nr(i) = data(i).nr; nc(i) = data(i).nc;
        sec(i) = data(i).sec; ser(i) = data(i).ser;
        xr(i) = data(i).xr; xc(i) = data(i).xc;
        
        if isempty(nc(i)); nc(i) = nan(); end
        
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


function [r2 rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of
% R-square, based on comparing the variability of the estimation errors
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
%
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true;
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1};
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
        % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));
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


function ciplot(lower,upper,x,colour);
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<4
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

fill([x fliplr(x)],[upper fliplr(lower)],colour, 'LineStyle', 'none')
end


function [beta, keepInd] = rmvOutlier (beta, keepInd);
% remove values that deviate from median by more than 3 median absolute
% deviations
try; keepInd; catch; keepInd = ones(length(beta),1); end


MAD = median(abs(beta - median(beta)));
keepInd((abs(beta) - abs(median(beta))) > 3*MAD) = 0;
if sum(keepInd) ~= 0
    beta = beta(keepInd == 1);
end

end


