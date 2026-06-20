
function [yrecon, stats] = nested_Reconstruction(Y, X, cal, ver, method, varargin)
% Combined reconstruction approach
% Inputs
% Y: Climate time-series
% X: Proxy array (same number of rows as Y), with each column representing a proxy series
% cal: Index of calibration set
% ver: Index of verification set
% method: Reconstruction method ('CPS', 'PCR', 'PCA')
% weightType: Type of weight to be used for CPS method ('r', 'r2', 'rp', or 'none')
% Outputs
% yrecon: Proxy-reconstructed climate time-series
% stats: Structure with reconstruction statistics

    [n, p] = size(X);
    Y_cal = nan(n,1);
    Y_cal(cal) = Y(cal);

    yrecon = nan(size(Y_cal));

    stats = struct('cal_r', nan(n,1), 'cal_rmse', nan(n,1), ...
                   'ver_r', nan(n,1), 'ver_rmse', nan(n,1), ...
                   'ver_CE', nan(n,1), 'ver_RE', nan(n,1));
    fieldsToUpdate = {'ver_r', 'ver_rmse', 'ver_CE', 'ver_RE'};

    % Sort proxies based on end times
    endTimes = arrayfun(@(col) find(~isnan(X(:, col)), 1, 'last'), 1:p);
    [~, sortOrder] = sort(endTimes, 'ascend'); % Sort from closest to farthest end time
    sortedX = X(:, sortOrder);

    % nesting approach
    for i = 1:p
        if i == 1
            % Select the reconstruction method based on the input 'method'
            [ yrecon, statsTemp ] = performReconstruction(Y, sortedX, cal, ver, method, varargin{:});

            a = find(all(~isnan(sortedX(:,:)), 2) == 1, 1, 'last' );
            yrecon(a+1:end) = NaN;
            stats = updateStats(stats, 1, a, statsTemp); % 1:a
            
            if all(isnan(yrecon(ver(1:a))))
                for fieldName = fieldsToUpdate
                    stats.(fieldName{1})(1:a) = NaN;
                end
            end

        else
            % Select the reconstruction method based on the input 'method'
            if i == p
                [currentRecon, statsTemp ] = performReconstruction(Y, sortedX(:,i), cal, ver, method, varargin{:});
            else
                [currentRecon, statsTemp ] = performReconstruction(Y, sortedX(:,i:end), cal, ver, method, varargin{:});
            end

            % common_period = ~isnan(yrecon) & ~isnan(sortedX(:,i-1)); % find time range of X(:,i-1)
            % currentRecon = scale(Y_calib, currentRecon, common_period);

            % Combine the current reconstruction with the final reconstruction
            a = find(~isnan(sortedX(:,i-1)) == 1, 1, 'last' );
            b = find(~isnan(sortedX(:,i)) == 1, 1, 'last' );

            yrecon(a+1:b) = currentRecon(a+1:b);
            stats = updateStats(stats, a+1, b, statsTemp); % a+1:b
            
            c = find(all(~isnan(sortedX(:,i:end)),2) == 1, 1 );
            if all(isnan(yrecon(ver(c:b))))
                for fieldName = fieldsToUpdate
                    stats.(fieldName{1})(c:b) = NaN;
                end
            end

        end
    end

    % Sort proxies based on start times
    startTimes = arrayfun(@(col) find(~isnan(X(:, col)), 1, 'first'), 1:p);
    [~, sortOrder] = sort(startTimes, 'descend'); % Sort from latest to closest start time
    sortedX = X(:, sortOrder);
    
    % nesting approach
    for i = 1:p
        if i == 1
            a = find(all(~isnan(sortedX(:,:)), 2) == 1, 1 );
            yrecon(1:a-1) = NaN;
            fieldNames = fieldnames(stats);
            for j = 1:length(fieldNames)
                stats.(fieldNames{j})(1:a-1) = NaN;
            end
        else
            % Select the reconstruction method based on the input 'method'
            if i == p
                [currentRecon, statsTemp ] = performReconstruction(Y, sortedX(:,i), cal, ver, method, varargin{:});
            else
                [currentRecon, statsTemp ] = performReconstruction(Y, sortedX(:,i:end), cal, ver, method, varargin{:});
            end

            % common_period = ~isnan(yrecon) & ~isnan(sortedX(:,i-1)); % find time range of X(:,i-1)
            % currentRecon = scale(Y_calib, currentRecon, common_period);

            % Combine the current reconstruction with the final reconstruction
            a = find(~isnan(sortedX(:,i)) == 1, 1 );
            b = find(~isnan(sortedX(:,i-1)) == 1, 1 )-1;
            yrecon(a:b) = currentRecon(a:b);
            stats = updateStats(stats, a, b, statsTemp); % a+1:b

            c = find(all(~isnan(sortedX(:,i:end)),2) == 1, 1, 'last' );
            if all(isnan(yrecon(ver(a:c))))
                for fieldName = fieldsToUpdate
                    stats.(fieldName{1})(a:c) = NaN;
                end
            end

        end
        
        % Display progress
        fprintf('Completed reconstruction for proxy %d/%d\n', i, p);
    end

    % Reconstruct the climate time-series by scaling with target's mean and standard deviation
    ind = ~isnan(yrecon) & ~isnan(Y_cal);
    yrecon = scale(Y_cal,yrecon,ind);

end

function stats = updateStats(stats, a, b, statsTemp)
        % calibration
        stats.cal_r(a:b) = statsTemp.cal_r;
        stats.cal_rmse(a:b) = statsTemp.cal_rmse;
        % verification
        stats.ver_r(a:b) = statsTemp.ver_r;
        stats.ver_rmse(a:b) = statsTemp.ver_rmse;
        stats.ver_CE(a:b) = statsTemp.ver_CE;
        stats.ver_RE(a:b) = statsTemp.ver_RE;
end

% Select the reconstruction method based on the input 'method'
function [ Recon, stats ] = performReconstruction(Y, X, cal, ver, method, varargin)
    % Check if weightType is provided for methods that require it
    if strcmp(method, 'CPS') && isempty(varargin)
        error('WeightType argument is required for CPS method.');
    end
    
    switch method
        case 'CPS'
            weightType = varargin{1};
            [ Recon, stats ] = cps_weighted(Y, X, cal, ver, weightType);
        case 'PCA'
            [ Recon, stats ] = pca_Reconstruction(Y, X, cal, ver);
        case 'PCR'
            [ Recon, stats ] = pcr_Reconstruction(Y, X, cal, ver);
        case 'RegEM'
            % normalized
            %%%%%%%%%%%%
            options = varargin{1}; % Options for RegEM
            [ Recon, stats] = regem_perform(Y, X, cal, ver, options);
        otherwise
            error('Invalid reconstruction method. Choose ''CPS'', ''PCR'', or ''PCA''.');
    end
end