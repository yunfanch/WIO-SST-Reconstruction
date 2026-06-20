function [ yrecon, stats ] = cps_weighted( Y, X, cal, ver, weightType )
% Composite-plus-scale: proxies weighted with scaling
% Inputs
% Y_calib: Climate time-series in calibration period, filled with NaNs for years to reconstruct
% X: Proxy array (same number of rows as y), with each column representing a proxy series
% weightType: Type of weight to be used ('r', 'r2', 'rp', or 'none')
% Outputs
% yrecon: proxy-reconstructed climate time-series
% stats: structure with Pearson correlation, p-value, variance-scaled R2 (McCarroll et al. 2015), and RMSE

    Y_cal = nan(size(Y));
    Y_cal(cal) = Y(cal);

    if nargin < 3
        % Default to 'rp' (r*(1-p)) if weightType is not provided
        weightType = 'rp';
    end
    
    [~, n] = size(X);
    
    % Calculate Pearson correlation and p-values
    [r, p] = corr(X, Y_cal, 'rows','pairwise');
    
    % Find indices of columns with positive correlation (r)
    positiveRIndices = find(r > 0);
    
    % Print the indices in the command window
    % fprintf('Columns with positive correlation (r > 0): %s\n', num2str(positiveRIndices));
    
    % Choose weights based on weightType
    switch weightType
        case 'r'
            weights = r';
        case 'r2'
            weights = r'.^2.*sign(r');
        case 'rp'
            weights = r'.*(1-p');
        case 'none'
            weights = ones(1, n).*sign(r'); % No weights
        otherwise
            error('Invalid weightType. Use ''r'', ''r2'', ''rp'', or ''none''.');
    end
    
    % % Select only columns with negative correlation
    % X = X(:, r <= 0);
    % weights = weights(r <= 0);
    
    if any(r < 0)
        % Normalize each proxy array 
        X_normalized = bsxfun(@rdivide, bsxfun(@minus, X, mean(X,'omitmissing')), std(X,'omitmissing'));
        
        % Apply weights to the normalized proxy array
        weighted_X = bsxfun(@times, X_normalized, weights);
        
        % Calculate the weighted composite index
        idx = sum(weighted_X, 2,'omitnan');
        
        % Reconstruct the climate time-series by scaling with target's mean and standard deviation
        ind = ~isnan(idx) & ~isnan(Y_cal);
        yrecon = scale(Y_cal,idx,ind);
    else
        yrecon = Y_cal;
    end

    % Statistical Tests Used to Assess the Coral Reconstruction
    % using an independent calibration-validation tests 
    [stats.cal_r, stats.cal_p] = corr(Y(cal), yrecon(cal), 'rows','pairwise');
    stats.cal_r2vs = 2 * abs(stats.cal_r) - 1;
    stats.cal_rmse = sqrt(mean((Y(cal) - yrecon(cal)).^2, 'omitnan'));
    stats.cal_mae = mean(abs(Y(cal) - yrecon(cal)), 'omitnan');

    [stats.ver_r, stats.ver_p] = corr(Y(ver), yrecon(ver), 'rows','pairwise');
    stats.ver_rmse = sqrt(mean((Y(ver) - yrecon(ver)).^2, 'omitnan'));
    stats.ver_CE = CE( Y, yrecon, ver);
    stats.ver_RE = RE( Y, yrecon, cal, ver);

end
