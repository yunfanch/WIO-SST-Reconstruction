function [y, t] = aggregate(x, time, method)
    % Validate inputs
    if ~isa(time, 'datetime')
        error('Time input must be of datetime type.');
    end
    if numel(x) ~= numel(time)
        error('Data and time inputs must have the same length.');
    end
    
    % Initialize outputs
    y = [];
    t = datetime([], [], []);

    % Extract year, month for further processing
    [yr, mo] = ymd(time);
    
    switch lower(method)
        case 'bimonthly'
            % Group by bimonthly periods
            bimonthIdx = ceil(mo/2) + (yr - min(yr))*12/2;
            
        case 'seasonal'
            % Define seasons: DJF(12,1,2) - index1, MAM(3,4,5) - index2,
            % JJA(6,7,8) - index3, SON(9,10,11) - index4
            % Mapping months to their respective season groups
            seasonIdx = ceil(mod(mo+1,12)/3);
            seasonIdx(seasonIdx == 0) = 4;
            seasonIdx(mo == 12) = seasonIdx(mo == 12) + 6;
            seasonIdx = seasonIdx + (yr - min(yr))*12/2;

        case 'annual'
            annualIdx = yr;
       
        otherwise
            error('Invalid method. Choose between "bimonthly", "season", or "annual".');
    end

    % Handle different aggregation methods
    if strcmpi(method, 'annual')
        uniqueYears = unique(annualIdx);
        for i = 1:length(uniqueYears)
            idx = annualIdx == uniqueYears(i);
            y(end+1) = mean(x(idx));
            t(end+1) = datetime(uniqueYears(i), 7, 1); % Use July 1st as representative date
        end
    elseif strcmpi(method, 'seasonal')
        uniqueSeasons = unique(seasonIdx);
        for i = 1:length(uniqueSeasons)
            idx = seasonIdx == uniqueSeasons(i);
            y(end+1) = mean(x(idx));
            % Find the median date for the period
            t(end+1) = time(find(idx,1,'first')) + days(floor((median(days(time(idx) - time(find(idx,1,'first')))))));
        end
    else strcmpi(method, 'biomonthly')
        uniqueMonths = unique(bimonthIdx);
        for i = 1:length(uniqueMonths)
            idx = bimonthIdx == uniqueMonths(i);
            y(end+1) = mean(x(idx));
            % Find the median date for the period
            t(end+1) = time(find(idx,1,'first')) + days(floor((median(days(time(idx) - time(find(idx,1,'first')))))));
        end
    end
end
