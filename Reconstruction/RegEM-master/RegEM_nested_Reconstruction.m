function [yrecon, stats ] = RegEM_nested_Reconstruction(Y, X, cal, ver, options)
% the ‘nesting’ reconstruction approach
% Inputs
% Y: Climate time-series
% X: Proxy array (same number of rows as Y), with each column representing a proxy series
% cal: Index of calibration set
% ver: Index of verification set
% Outputs
% Recon: proxy-reconstructed climate time-series
% stats: structure with Pearson correlation, p-value, variance-scaled R2 (McCarroll et al. 2015), and RMSE

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
            % apply RegEM
            [ yrecon, statsTemp ] = regem_cfr( Y, sortedX, cal, ver, options);

            a = find(all(~isnan(sortedX(:,:)), 2) == 1, 1, 'last' );
            yrecon(a+1:end) = NaN;
            stats = updateStats(stats, 1, a, statsTemp); % 1:a

            if all(isnan(yrecon(ver(1:a))))
                for fieldName = fieldsToUpdate
                    stats.(fieldName{1})(1:a) = NaN;
                end
            end

        else
            % apply RegEM
            if i == p
                [ currentRecon, statsTemp] = regem_cfr( Y, sortedX(:,i), cal, ver, options);
            else
                [ currentRecon, statsTemp] = regem_cfr( Y, sortedX(:,i:end), cal, ver, options);
            end

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
            % apply RegEM
            if i == p
                [currentRecon, statsTemp ] = regem_cfr(Y, sortedX(:,i), cal, ver, options);
            else
                [currentRecon, statsTemp ] = regem_cfr(Y, sortedX(:,i:end), cal, ver, options);
            end
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