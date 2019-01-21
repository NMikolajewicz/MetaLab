function [output] = prepData_GUI_170917(input)
clearvars('-except', 'input');

% PREPDATA Import meta-analysis data from Microsoft Excel spreadsheet file

%   [DATACOLLECTION] = PREPDATA(FILE, SHEET, COLLECTIONNAME, COLLAPSE, SAVESET, MINGROUPSIZE)
%   reads data from FILE from specified worksheet SHEET and
%   returns meta-analytic input data as data structure DATACOLLETION.
%   COLLECTIONNAME specifies name of DATACOLLECTION. Optionally, COLLAPSE
%   specifies whether data is pooled according. Optionally, data
%   collection is saved if SAVESET is true.

%   Input Arguments:
%------------------------
%
%   FILE    String that specifies the name of the file to read. (ex. 'filename.xlsx')
%   SHEET   String that contains the worksheet name. (ex. 'Sheet1')
%           *Excel sheet must, at minimum, contain one of the following set of headers:
%           (i) 'ID', 'Study', 'xr', 'ser', 'nr'
%           (ii) 'ID, 'Study', 'xr', 'ser', 'nr', 'xc', 'sec', 'nc'
%           *Additional headers are optional, and are treated as covariates
%   COLLETIONNAME
%           String that specifies name of output data collection (ex. 'myDataCollection')
%   COLLAPSE
%           Logical (TRUE/FALSE). Specfies whether data is pooled.
%           *Pooled data indices must be specified numerically in worksheet
%           under header 'collapse'.
%           COLLETIONNAME
%   MINGROUPSIZE
%           Specifies minimum subgroup size.
%           *MINGROUPSIZE must be positive integer
%
%   last update: 14.09.17

try;
try; properties.file = input.file; catch; error('input file not specified'); end
try; properties.sheet = input.sheet{1}; catch; error('sheet name not specified'); end
try; properties.collectionName = input.collectionName{1}; catch; end
try; properties.collapse =input.collapse; catch; properties.collapse = false; end
try; properties.saveSet = input.saveSet; catch;properties.saveSet = false; end
try; properties.minGroupSize = input.minGroupSize; catch; properties.minGroupSize = 2; end

[num,txt,raw]=xlsread(properties.file, properties.sheet);
inputHeaders = raw(1,:);

reqHeaders = { 'ID', 'xr', 'ser', 'nr'};
reqHeaders_basal = {'xc', 'sec', 'nc'};


%% specify data of interest
n = 1;
for i = 1:length(reqHeaders)                            %assign data to structure
    for j = 1:length(inputHeaders)
        if strcmp(reqHeaders{i}, inputHeaders{j});
            necInd(n) = j; n = n + 1;
            for k = 2:size(raw,1)
                data(k-1).(reqHeaders{i}) = raw{k,j};
            end
        end
    end
end

for i = 1:length(reqHeaders_basal)                            %assign data to structure
    for j = 1:length(inputHeaders)
        if strcmp(reqHeaders_basal{i}, inputHeaders{j});
            necInd(n) = j; n = n + 1;
            for k = 2:size(raw,1)
                data(k-1).(reqHeaders_basal{i}) = raw{k,j};
            end
        end
    end
end

% check if collapse array is present

for j = 1:length(inputHeaders)
    if strcmp('collapse', inputHeaders{j});
        necInd(n) = j; n = n + 1;
        for k = 2:size(raw,1)
            data(k-1).collapse = raw{k,j};
        end
    end
end

% remove studyNames so they are not included as covariates
for j = 1:length(inputHeaders)
    if strcmp('Study', inputHeaders{j});
        necInd(n) = j; n = n + 1;
    end
end

% find available covariates
n = 1;
for i = 1:length(inputHeaders)
    if ~any(i == necInd)
        covariates{n} = inputHeaders{i};
        n = n+1;
    end
end

% add covariates to data structure
try;
for i = 1:length(covariates)
    for j = 1:length(inputHeaders)
        if strcmp(covariates{i}, inputHeaders{j});
            for k = 2:size(raw,1)
                data(k-1).( covariates{i}(find(~isspace(covariates{i})))) = raw{k,j};
            end
        end
    end
    covariates{i} =  covariates{i}(find(~isspace(covariates{i})));
end
catch; display('no covariates found'); covariates{1} = []; end

% remove nan, empty or non numeric entries
n = 1;
for i = 1:length(data)
    if ~isnan(data(i).xr) & ~isnan(data(i).ser) & isnumeric(data(i).xr) & isnumeric(data(i).ser) &  ~isempty(data(i).xr) & ~isempty(data(i).ser) & data(i).ser ~= 0
        testSet(n) = data(i); n = n+1;
    end
end

%% collapse studies with multiple data sets into single
if properties.collapse
    ident = [testSet.collapse]; % find all
    [uniqID,ia, ic] = unique(ident);
    
    n = 1;
    for i = 1:length(ia);
        xr{i} = []; ser{i} = []; Nr{i} = []; xc{i} = []; sec{i} = []; Nc{i} = [];
        for co = 1:length(covariates)
            cova(i).(covariates{co}) = [];
        end
        for j = 1:length(testSet);
            if testSet(j).collapse == ident(ia(i));
                xr{i} = [xr{i} testSet(j).xr];
                ser{i} = [ser{i} testSet(j).ser];
                Nr{i} = [ Nr{i} testSet(j).nr];
                try;
                    xc{i} = [xc{i} testSet(j).xc];
                    sec{i} = [sec{i} testSet(j).sec];
                    Nc{i} = [ Nc{i} testSet(j).nc];
                catch;
                end
                try; ind(i) = [ind(i) j]; catch; ind(i) =  j ; end
                for co = 1:length(covariates)
                    cova(i).(covariates{co}) = [cova(i).(covariates{co}) testSet(j).(covariates{co})];
                end
            end
        end
        
        %% pool intrastudy data sets together
        xp(i) = sum( xr{i} .* (Nr{i})) / sum(Nr{i});
        try;  xcp(i) = sum( xc{i} .* (Nc{i})) / sum(Nc{i}); catch; end
        vp(i) = sum((Nr{i}-1) .* (ser{i}.^2)) / sum((Nr{i}-1));
        try; vcp(i) = sum((Nc{i}-1) .* (sec{i}.^2)) / sum((Nc{i}-1)); catch; end
        Np(i) = sum(Nr{i});
        try;  Ncp(i) = sum(Nc{i}); catch; end;
        
        %% pool covariate data
        for jo = 1:length(covariates)
            if length(unique(cova(i).(covariates{jo}))) == 1
                finalCo(i).(covariates{jo}) = unique(cova(i).(covariates{jo}));
            else
                finalCo(i).(covariates{jo}) = [];
            end
        end
        testSet_uniq(i) = testSet(ind(i));
    end
    % consolidate FINAL_uniq set with collapsed data sets
    n = 1;
    for i = 1:length(testSet_uniq);
        FINAL(i).ID = uniqID(i);
        FINAL(i).xr = xp(i);
        FINAL(i).ser = sqrt(vp(i));
        FINAL(i).nr = Np(i);
        FINAL(i).xc = xcp(i);
        FINAL(i).sec = sqrt(vcp(i));
        FINAL(i).nc = Ncp(i);
        for jo = 1:length(covariates)
            FINAL(i).(covariates{jo})= finalCo(i).(covariates{jo});
        end
    end
    testSet = []; testSet = FINAL;
end


%% stratify & prepare data
%--------------------------------------------------------------------------
origSet = testSet;
stratChoice = menu('Subgroup Stratification (number of subgroup stratifications)', '1 level', '2 levels', '3 levels', '4 levels', '5 levels', 'no stratification (recommended)');
if stratChoice ~= 6
    for p = 1:stratChoice
        choice(p) = menu('define data set',covariates{:});
        display(['level ' num2str(p) ': ' covariates{choice(p)}]);
    end
    for p = 1:length(choice)
        subgroupCode = [];
        foc{p} = covariates{choice(p)};
        subgroupCode = [testSet.(foc{p})];
        uniqCode{p} = unique(subgroupCode);
        [uniqCode{p}, ~] = prepareCurveData(uniqCode{p}, uniqCode{p});
    end
    if length(choice) ==1;
        combinationVec = uniqCode{1};
    elseif length(choice) == 2;
        combinationVec = allcomb(uniqCode{1}, uniqCode{2});
    elseif length(choice) == 3
        combinationVec = allcomb(uniqCode{1}, uniqCode{2},uniqCode{3} );
    elseif length(choice) == 4
        combinationVec = allcomb(uniqCode{1}, uniqCode{2},uniqCode{3},uniqCode{4});
    elseif length(choice) == 5
        combinationVec = allcomb(uniqCode{1}, uniqCode{2},uniqCode{3},uniqCode{4},uniqCode{5} );
    end
else
    combinationVec = nan();
end
assignin('base', 'combinationVec', combinationVec);
if ~isempty(combinationVec)
    for nDS = 1:size(combinationVec, 1)
        clear ('testSet');  testSet = origSet;
        if ~isnan(combinationVec) % if covariates were selected for analysis
            for i = 1:length(testSet)
                incomplete(i) = 0;
                if (testSet(i).(foc{1}) ~= combinationVec(nDS,1)); incomplete(i) = 1; end;
                if length(choice) > 1
                    if (testSet(i).(foc{2}) ~= combinationVec(nDS,2)); incomplete(i) = 1; end; end;
                if length(choice) > 2
                    if (testSet(i).(foc{3}) ~= combinationVec(nDS,3)); incomplete(i) = 1; end; end;
                if length(choice) > 3
                    if (testSet(i).(foc{4}) ~= combinationVec(nDS,4)); incomplete(i) = 1; end; end;
                if length(choice) >4
                    if (testSet(i).(foc{5}) ~= combinationVec(nDS,5)); incomplete(i) = 1; end; end;
                testSet(i).incomplete = incomplete(i);
            end
            setName{1} = [foc{1} num2str(combinationVec(nDS,1))];
            for p = 2:length(choice)
                setName{1} = [setName{1} '-' foc{p}  num2str(combinationVec(nDS,p))];
            end
        else
            for i = 1:length(testSet) % if total set was selected for analysis
                incomplete(i) = 0;
                testSet(i).incomplete = incomplete(i);
            end
            setName{1} = 'totalSet';
        end
        %--------------------------------------------------------------------------
        
        n = 1;
        temp = testSet; clear ('testSet');
        
        for i = 1:length(temp);
            if temp(i).incomplete == 0
                testSet(n) = temp(i);
                n = 1 + n;
            end
        end
        
        try;
            FINAL_uniq = testSet;
        catch; FINAL_uniq = [];
        end
        
        try; FINAL_uniq = rmfield(FINAL_uniq, 'incomplete'); catch; end
        try; FINAL_uniq = rmfield(FINAL_uniq, 'collapse'); catch; end
        if ~isempty(FINAL_uniq)
            if length(FINAL_uniq) >= properties.minGroupSize
                try; load(properties.collectionName); catch; display('new data collection created'); end
                try; current = length(D)+1; catch; current = 1; end
                D(current).description = setName;
                D(current).data = FINAL_uniq;
                D(current).covariates = covariates;
                if properties.saveSet
                    save(properties.collectionName, 'D');
                    clear ('FINALSET'); clear ('testSet'); clear ('testSet_uniq'); clear ('nameHolder');
                end
            end
        end
    end
    
    
    
    assignin('base','FINAL_uniq',FINAL_uniq);
    try; setNames{1} = [D(1).description{1}]; catch; display('no data set stored, problem with input data'); end
    for i = 1:length(D)-1;
        setNames{1} = [setNames{1} ', ' (D(i+1).description{1})];
    end
    
    display(' ');
    display('Data preparation successful.')
    display(['Current data collection "' properties.collectionName '" contains '  num2str(length(D)) ' datasets']);
    display(['Stored data set: ' setNames{1}]);
    display(' ');
    display(['Use "' properties.collectionName '" as input for HETEROGENEITY, METAANALYSIS, or METAREGRESSION modules']);
    display(' ');
    display('Tip: run "prepData" again to add more datasets to current data collection');
    
    output = D;
    
    
    msgbox({['Current data collection "' properties.collectionName '" contains '  num2str(length(D)) ' datasets'],...
        ['Stored data set: ' setNames{1}],...
        ' ',...
        ['Use "' properties.collectionName '" as input for HETEROGENEITY, METAANALYSIS, or METAREGRESSION modules'],...
        ' ',...
        'Tip: run "prepData" again to add more datasets to current data collection'}, 'Data preparation successful');
    clearvars('-except', 'output', 'covariates', 'checkThis');
    
else
    try;
        try; load(properties.collectionName); catch; end
        display(' ');
        display('There was insuccifient data avaiable for further analysis.')
        display(['Current data collection "' properties.collectionName '" contains '  num2str(length(D)) ' datasets']);
        display(['Stored data set: ' setNames{1}]);
        display(' ');
        display(['Use "' properties.collectionName '" as input for HETEROGENEITY, METAANALYSIS, or METAREGRESSION modules']);
        display(' ');
        display('Tip: run "prepData" again to add more datasets to current data collection');
        
        output = D;
                
        msgbox({['Current data collection "' properties.collectionName '" contains '  num2str(length(D)) ' datasets'],...
            ['Stored data set: ' setNames{1}],...
            ' ',...
            ['Use "' properties.collectionName '" as input for HETEROGENEITY, METAANALYSIS, or METAREGRESSION modules'],...
            ' ',...
            'Tip: run "prepData" again to add more datasets to current data collection'}, 'Insuccifient data available for further analysis');
    catch;
        msgbox({'There was insuccifient data avaiable for further analysis.'});
        output = D;
    end
end

catch;     
  msgbox({'Error while running prepare data module',...
      '                                                              ',...
      'Ensure inputs are complete and correct'},'Error', 'Error');  
  output = [];
end
end
function A = allcomb(varargin)

% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
%    of the number of elements of the N inputs. This functionality is also
%    known as the Cartesian Product. The arguments can be numerical and/or
%    characters, or they can be cell arrays.

% Tested in Matlab R2015a
% version 4.1 (feb 2016)
% (c) Jos van der Geest
% email: samelinoa@gmail.com


narginchk(1,Inf) ;

NC = nargin ;

% check if we should flip the order
if ischar(varargin{end}) && (strcmpi(varargin{end},'matlab') || strcmpi(varargin{end},'john')),
    % based on a suggestion by JD on the FEX
    NC = NC-1 ;
    ii = 1:NC ; % now first argument will change fastest
else
    % default: enter arguments backwards, so last one (AN) is changing fastest
    ii = NC:-1:1 ;
end

args = varargin(1:NC) ;
% check for empty inputs
if any(cellfun('isempty',args)),
    warning('ALLCOMB:EmptyInput','One of more empty inputs result in an empty output.') ;
    A = zeros(0,NC) ;
elseif NC > 1
    isCellInput = cellfun(@iscell,args) ;
    if any(isCellInput)
        if ~all(isCellInput)
            error('ALLCOMB:InvalidCellInput', ...
                'For cell input, all arguments should be cell arrays.') ;
        end
        % for cell input, we use to indices to get all combinations
        ix = cellfun(@(c) 1:numel(c), args,'un',0) ;
        
        % flip using ii if last column is changing fastest
        [ix{ii}] = ndgrid(ix{ii}) ;
        
        A = cell(numel(ix{1}),NC) ; % pre-allocate the output
        for k=1:NC,
            % combine
            A(:,k) = reshape(args{k}(ix{k}),[],1) ;
        end
    else
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}),[],NC) ;
    end
elseif NC==1,
    A = args{1}(:) ; % nothing to combine
    
else % NC==0, there was only the 'matlab' flag argument
    A = zeros(0,0) ; % nothing
end

end

