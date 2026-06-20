function [ yrecon, stats ] = pcr_Reconstruction(Y, X, cal, ver)
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
    %[n,numPCs] = size(X_normalized);

    % Perform PCA
    common_period = all(~isnan([Y_cal X_normalized]), 2);
    xx = X_normalized(common_period,:);
    yy = Y_cal(common_period);
    [n,m] = size(xx);
    [PCALoadings, PCAScores, D] = pca(xx);

    % PCALoadings: This represents the principal component loadings or eigenvectors. 
    % Each column of PCALoadings corresponds to a principal component, 
    % and it shows the weights (loadings) of the original variables on each principal component. 
    % These loadings help you understand the contribution of each original variable to the principal components.
    % PCAScores: These are the scores or projections of your data onto the principal components. 
    % Each row of PCAScores corresponds to an observation (data point), 
    % and each column corresponds to a principal component. 
    % PCAScores allows you to represent your original data in the principal component space.
    % D: This represents the eigenvalues of the covariance matrix of your data. 
    % The eigenvalues are a measure of the variance explained by each principal component. 
    % They are useful for assessing the importance or contribution of each principal component to the total variance in your data.

    % Find the optimal number of principal components
    % crossval combined with a simple function to compute the sum of squared errors for PCR, 
    % can estimate the MSPE,  using 10-fold cross-validation.
    PCRmsep = sum(crossval(@pcrsse,xx,yy,'KFold',10),1) / n;
    %figure;plot(0:7,PCRmsep,'r-^');
    [~, optimalPCs] = min(PCRmsep);
    optimalPCs = optimalPCs-1; % PCRmsep from 0 to numPCs
    % optimalPCs = m;
    % figure; plot(0:m,PCRmsep,'r-^');

    % Select the optimal number of principal components
    betaPCR = regress(yy-mean(yy), PCAScores(:,1:optimalPCs));
    betaPCR = PCALoadings(:,1:optimalPCs)*betaPCR;
    betaPCR = [mean(yy) - mean(xx)*betaPCR; betaPCR];
    yPredicted = [ones(size(X_normalized,1),1) X_normalized]*betaPCR;
    % figure;plot(Y_calib,Y_calib,'k',Y_calib,yPredicted,'r^');
             
    % Reconstruct the climate time-series by scaling with target's mean and standard deviation
    ind = ~isnan(yPredicted) & ~isnan(Y_cal);
    yrecon = scale(Y_cal,yPredicted,ind);

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
