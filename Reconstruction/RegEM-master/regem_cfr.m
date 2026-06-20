function [field_r,diagn] = regem_cfr(field,proxy,calib,options)
% Function [field_r,diagn] = regem_cfr(field,proxy,calib,options)
%  Performs RegEM-based climate reconstruction 
%
%  inputs:  - field, 2D climate field (nf x pf)
%           - proxy, proxy matrix (np x pp)
%           - calib, index of calibration set (< np)
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
[nf,pf] = size(field); % field dimensions
[np,pp] = size(proxy); % proxy matrix dimensions

%	Assemble climate field/proxy data matrices
X_in  = NaN(np,pf+pp);
inst    = field(calib,:); % instrumental field
ni      = size(inst,1);
X_in(calib,1:pf)   = inst;   % Put in instrumental data
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
    X_hat = X(calib, avail{j}) * B{j};
    Xp    = X_hat(:,1:pf);
else
    warning('Prediction over the calibration period could not be performed')
    warning('Check that proxy availability does not vary at t = min(ti)')
    Xp    = X(:,1:pf);
end

% add mean to centered data matrices
Xp = Xp + repmat(M(1:pf), [ni 1]);
% assign predicted values of the field over the calibration interval
field_r(calib,:) = Xp(:,1:pf);

% give peff a more user-friendly form
peff_e = zeros(np,1);
for k = 1:npat
   if ~isempty(peff{k})
      peff_e(iptrn==k) = mean(peff{k});
	else
		peff_e(iptrn==k) = peff{k}; % assign peff from the penultimate pattern of missing values
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

if insample
   [Xp,RE,R2] = insample_pred_regem(X,M,B,[1:pf],[pf+1:pp],calib,avail,iptrn);
   diagn.RE = RE;
   diagn.R2 = R2;  
end


end
