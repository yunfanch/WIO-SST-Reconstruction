function result = regem_Reconstruction(X, options, tracker)
%REGEM Multiproxy climate reconstruction with Regularized 
%   Expectation-Maximization
%   result = regem(data, options, tracker) calculates the climate signal
%   that best matches the pairwise comparisons of all proxy data. 
%
%   If options is omitted, defaults will be used. If tracker is omitted, no
%   output is generated.
%   
%   If data is set to 'info', result returns default options and
%   information about their meaning.

if nargin == 0
    error('RECON:REGEM','At least one input parameter has to be defined.');    
end

% Check options and fill in values that are not defined
if (~exist('options','var'))
    options = [];
end
options = defaults(options);  

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

% Current implementation is optimized for TTLS and no other regression
% method is supported 
options.regress = 'ttls';

%% Preprocess

[n, p]       = size(X);
% number of degrees of freedom for estimation of covariance matrix
dofC         = n - 1;            % use degrees of freedom correction      
% if X is a vector, make sure it is a column vector (a single variable)

% get indices of missing values and initialize matrix of imputed values
missing = isnan(X);
nrmissing = nnz(missing);
if nrmissing == 0
    warning('No missing value flags found.')
    return                                      % no missing values
end
[~,kmis]  = ind2sub([n, p], find(missing));

% default values
Xmis = nan(n,p);
Xerr = inf(n,p);
options.stagTolerance = 1e-6;
options.maxIterations = 1e4;
options.neigs  = min(n-1, p);
options.inflation    = 1;

% initial estimates of missing values
[X, M] = center(X);
X(missing)  = zeros(nrmissing, 1);   % fill missing entries with zeros

% initial estimate of covariance matrix      
C = X'*X / dofC;      

  
if strcmpi(options.regress,'ttls')
   [~,S,~] = svd(X./repmat(diag(C)',size(X,1),1));
   trunc = eigenselect(diag(S).^2);
end

% Find patterns of missing data
[missingPatterns,~,rowToMisPattern] = unique(missing,'rows');
firstPattern = 1;
if (all(~missingPatterns(1,:)))
    % Full of data. Skip this in pttls
    firstPattern = 2;
end
iteration = 0;
rdXmis  = Inf;

tracker('RegEM');
while (iteration < options.maxIterations && rdXmis > options.stagTolerance)
    iteration = iteration + 1;

    % initialize for this iteration ...
    CovRes = zeros(p,p);       % ... residual covariance matrix
    peff_ave = 0;                % ... average effective number of variables 

    % Scale variables to unit variance
    dataStd = sqrt(diag(C)); 
    
    dataStd(abs(dataStd) < eps) = 1;   % do not scale constant variables
    X = X ./ repmat(dataStd',n,1);
    C = C ./ (dataStd*dataStd');    

    % compute eigendecomposition of correlation matrix
    [V, d]   = peigs(C, options.neigs);
    peff_ave = (dofC - trunc) * nrmissing;
    Sglobal = V(:,trunc+1:end);
    Sglobal = Sglobal * diag(d(trunc+1:end))* Sglobal';
    
    
    for patIndex=firstPattern:size(missingPatterns,1)         % cycle over records
        misMask = missingPatterns(patIndex,:);
        rowPattern = rowToMisPattern == patIndex;        
        % 'ttls'
        % truncated total least squares with fixed truncation parameter
        [B, S]   = pttls(V, d, misMask, trunc, Sglobal);

        dofS     = dofC - trunc;         % residual degrees of freedom

        % inflation of residual covariance matrix
        S        = options.inflation * S;

        % bias-corrected estimate of standard error in imputed values
        Xerr(rowPattern, misMask) = repmat(dofC/dofS * sqrt(diag(S))',nnz(rowPattern),1);
        

        % missing value estimates
        Xmis(rowPattern, misMask)   = X(rowPattern, ~misMask) * B;

        % add up contribution from residual covariance matrices
        inplaceadd(CovRes,S*nnz(rowPattern),misMask);
    end % loop over missing patterns
    
    % rescale variables to original scaling 
    D = repmat(dataStd',n,1);
    X  = X .* D;
    Xerr = Xerr .* D;
    Xmis = Xmis .* D;
%     C          = C .* repmat(D', p, 1) .* repmat(D, 1, p);
    CovRes = CovRes.*(dataStd*dataStd');

    % rms change of missing values
    dXmis = norm(Xmis(missing) - X(missing)) / sqrt(nrmissing);
    
    % relative change of missing values
    nXmis_pre  = norm(X(missing) + M(kmis)') / sqrt(nrmissing);    
    if nXmis_pre < eps
        rdXmis   = Inf;
    else
        rdXmis   = dXmis / nXmis_pre;
    end

    % update data matrix X
    X(missing)  = Xmis(missing);
    
    % re-center data and update mean
    [X, Mup]   = center(X);                  % re-center data
    M          = M + Mup;                    % updated mean vector

    % update covariance matrix estimate
    C          = (X'*X + CovRes)/dofC; 
    
%     figure(3); clf; plot(X(:,2:end)); hold all;
%     plot(X(:,1), 'linewidth',2,'color',[0.7 0.7 0.7]);
%     drawnow;
%     fprintf('%d: %5.9f\n',iteration, rdXmis);

    tracker('RegEM', iteration, options.maxIterations);
end                                        % EM iteration

% add mean to centered data matrix
X  = X + repmat(M, n, 1);  

result.field = X;
% result.field = X(:,1:size(data.instrumental.data,1))';


