function [ yrecon, stats ] = pca_Reconstruction(Y, X, cal, ver)
% Reconstructions based on PCA assume that the underlying gradient common to a group of time series 
% significantly correlated with climate time-series in the calibration window should be equivalent to climate time-series. 
% For PCA-based reconstructions, we therefore only used records that are correlated with climate time-series in the calibration window. 
% For PCA reconstructions, we performed PCA on the calibration window (that is, 
% in which we know that the proxy records are correlated with the climate time-series) 
% and then multiplied the loading of each proxy record on PC1 by the complete time series of the proxy records. 
% The contribution of each record to the PCA was weighted according to the strength of its correlation with climate time-series.
% The direction of a PC axis is arbitrary. To align the temporal subsets, we flipped (if necessary) PC1 of the calibration window 
% to make it positively correlated with the target index and then aligned PC1 of subsequent temporal subsets 
% to be positively correlated with their predecessors.       
% Inputs
% Y_calib: Climate time-series in calibration period, filled with NaNs for years to reconstruct
% X: Proxy array (same number of rows as y), with each column representing a proxy series
% Outputs
% yrecon: proxy-reconstructed climate time-series
% stats: structure with Pearson correlation, p-value, variance-scaled R2 (McCarroll et al. 2015), and RMSE

    Y_cal = nan(size(Y));
    Y_cal(cal) = Y(cal);

    % Normalize each proxy array 
    X_normalized = bsxfun(@rdivide, bsxfun(@minus, X, mean(X,'omitmissing')), std(X,'omitmissing'));
    
    % Perform PCA
    [PCALoadings, PCAScores, D] = pca(X_normalized);
          
    % PCALoadings (loadings matrix): Each column represents a principal component,
    % each row corresponds to a variable (or feature, proxy record) in coralProxies.
    % Each element in PCALoadings represents the weight or contribution of the original variable in the corresponding principal component.
    % The direction of a PC axis is arbitrary. 
    % Multiply the loading of each proxy record on PC1 by X_normalized
    pc1Loadings = PCALoadings(:, 1);
    PC_1 = X_normalized * pc1Loadings;     % Equivalent to the following calculation
    % Z=PCALoadings'*X_normalized';
    % PC_1=Z(1,:); 
    % or
    % the complete time series of the proxy records
    % weighted_X = X_normalized .* repmat(pc1Loadings', size(X_normalized, 1), 1);
    % Calculate the weighted composite index
    % PC_1 = sum (weighted_X,2,'omitnan');

    % To align the temporal subsets, we flipped (if necessary) PC1 of the calibration period 
    % to make it positively correlated with the target index 
    % and then aligned PC1 of subsequent temporal subsets to be positively correlated with their predecessors.
    [r, ~] = corr(PC_1, Y_cal, 'rows','pairwise');
    if r<0
        PC_1 = PC_1*(-1);
    end
            
    % Reconstruct the climate time-series by scaling with target's mean and standard deviation
    ind = ~isnan(PC_1) & ~isnan(Y_cal);
    yrecon = scale(Y_cal,PC_1,ind);
    
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
