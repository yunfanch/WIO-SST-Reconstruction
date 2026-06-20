function [yrecon,stats] = regem_perform(Y,X,cal,ver,options)
% Function [field_r,diagn] = regem_cfr(field,proxy,calib,options)
%  Performs RegEM-based climate reconstruction 
%
%  inputs:  - Y, 2D climate field (nf x pf)
%           - X, proxy matrix (np x pp)
%           - cal, index of calibration set (< np)
%           - ver, index of verification set (< np)
%           - options, structure of optional RegEM parameters. 
%
%  outputs: - field_r, reconstructed field (np x pf)
%           - diagn, structure of diagnostic outputs, including
%              * Xerr, estimate of imputation error 
%              * avail, vector of available values
%              * miss, vector of missing values
%              * iptrn, index of pattern to which each row belongs
%              * B, regression matrix
%              * RE,R2: "reduction of error" and R-squared in-sample statistics 
%                  (computed iff options.insample = 1).  
%              * peff, number of effective parameters
%
% All outputs are given for each pattern of missing values, except the last
% three: RE and R2 have dimensions (np x pf), and peff is a np x 1 vector.
%
% History: created  22-Nov-2013 09:07:55  by Julien Emile-Geay (USC)
%           edited  26-Nov-2013 13:57     by JEG to complete description and fix
%                definition of peff over the calibration interval
% ===========================================================================

% Process options
if nargin < 3 | isempty(options)
   fopts = [];
else
   fopts = fieldnames(options);
end

if strmatch('insample', fopts)  % compute in-sample skill or not?
   insample = options.insample;
else
   insample = false;
end

if strmatch('X0', fopts)  % initial guess for temperature
   X0_given = 1;
   X0 = options.X0;
else
   X0_given = 0;
end
%   ====== end options processsing  =======

% Define Time Parameters
% ========================
% normalized
for i = 1:size(X,2)
    X_mean(:,i) = mean(X(:,i),'omitmissing');
    X_std(:,i) = std(X(:,i),'omitmissing');
    X_nor(:,i) = (X(:,i)-X_mean(:,i))./X_std(:,i);
end
for i = 1:size(Y,2)
    Y_mean(:,i) = mean(Y(:,i),'omitmissing');
    Y_std(:,i) = std(Y(:,i),'omitmissing');
    Y_nor(:,i) = (Y(:,i)-Y_mean(:,i))./Y_std(:,i);
end

% Initialize Y_cal with NaN and assign values only for calibration indices
Y_cal = nan(size(Y_nor));
Y_cal(cal,:) = Y_nor(cal,:);

% Identify periods with available data in either Y_cal or X
recon_period = any(~isnan([Y_cal, X_nor]), 2);

% Extract relevant proxy data and corresponding field values
proxy = X_nor(recon_period,:)*(-1);
field = Y_cal(recon_period,:);
[nf,pf] = size(field); % field dimensions
[np,pp] = size(proxy); % proxy matrix dimensions

% Initialize the reconstruction output array
yrecon = nan(size(Y_nor));

% Determine the indices within the calibration period
cal_recon = false(size(Y_nor));
[~,ia,ib] = intersect(Y_cal,field);
cal_recon(ib) = true;

% Prepare the input matrix with calibration field values and proxy data
X_in  = NaN(np,pf+pp);
inst    = field(cal_recon,:); % instrumental field
ni      = size(inst,1);
X_in(cal_recon,1:pf)   = inst;   % Put in instrumental data
X_in(:,pf+1:pf+pp) = proxy;  % add (possibly incomplete) proxy data

%  Assign RegEM options
%==================================
if X0_given
%    X0 =   NaN(np,pf+pp);
%    X0(:,1:pf)    = field0;  % Put in temperature data
%    X0(:,pf+1:pf+pp) = proxy;  % the rest is (possibly incomplete) proxy data
   options.Xmis0 = X0;
end

% Apply RegEM 
[X, M, C, Xerr, B , peff, avail, miss, iptrn] = regem(X_in,options);

% extract reconstructed field
field_r = X(:,1:pf);

% compute predicted values of the field over the calibration interval
npat = length(avail); % number of patterns
Xp = zeros(ni,pf); %output array

% remove mean
X = X - repmat(M, [np 1]);

j = npat - 1; % no need to compute predictions using all patterns of missing values
if sum(isnan(B{j}))==0
    % make instrumental prediction
    X_hat = X(cal_recon, avail{j}) * B{j};
    Xp    = X_hat(:,1:pf);
else
    warning('Prediction over the calibration period could not be performed')
    warning('Check that proxy availability does not vary at t = min(ti)')
    Xp    = X(:,1:pf);
end

% add mean to centered data matrices
Xp = Xp + repmat(M(1:pf), [ni 1]);
% assign predicted values of the field over the calibration interval
field_r(cal_recon,:) = Xp(:,1:pf);

% give peff a more user-friendly form
peff_e = zeros(np,1);
for k = 1:npat
    if ~isempty(peff{k})
      peff_e(iptrn==k) = mean(peff{k});
	else
		peff_e(iptrn==k) = []; % assign peff from the penultimate pattern of missing values
	end	
end

% Assign output data structures
diagn.err = Xerr;
diagn.avail = avail;
diagn.miss  = miss;
diagn.peff  = peff_e;
diagn.iptrn = iptrn;
diagn.B     = B;
diagn.C     = C;
diagn.X     = X;
diagn.M     = M;

% if insample
%    [Xp,RE,R2] = insample_pred_regem(X,M,B,[1:pf],[pf+1:pp],cal,avail,iptrn);
%    diagn.RE = RE;
%    diagn.R2 = R2;  
% end


% Reconstruct the climate time-series by scaling with target's mean and standard deviation
yrecon(recon_period) = field_r;
yrecon = scale(Y,yrecon,cal);


    % Statistical Tests Used to Assess the Coral Reconstruction
    % using an independent calibration-validation tests 
    [stats.cal_r, stats.cal_p] = corr(Y(cal), yrecon(cal), 'rows','pairwise');
    stats.cal_r2vs = 2 * abs(stats.cal_r) - 1;
    stats.cal_rmse = sqrt(mean((Y(cal) - yrecon(cal)).^2, 'omitnan'));
    stats.cal_mae = mean(abs(Y(cal) - yrecon(cal)), 'omitnan');

    % fieldsToUpdate = {'ver_r', 'ver_rmse', 'ver_CE', 'ver_RE'};
    % if all(isnan(yrecon(ver)))
    %     for fieldName = fieldsToUpdate
    %         stats.(fieldName{1}) = NaN;
    %     end
    % else
        [stats.ver_r, stats.ver_p] = corr(Y(ver), yrecon(ver), 'rows','pairwise');
        stats.ver_rmse = sqrt(mean((Y(ver) - yrecon(ver)).^2, 'omitnan'));
        stats.ver_CE = CE( Y, yrecon, ver);
        stats.ver_RE = RE( Y, yrecon, cal, ver);
 %   end

end
