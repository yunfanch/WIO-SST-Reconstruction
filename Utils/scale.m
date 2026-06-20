function scaled_index = scale(reference, X, common_period)
% Scale the X to match the mean and variance of the reference_nest
% using common_period for alignment

    % Normalize the index
    mean_to_scale = mean(X(common_period),'omitmissing');
    std_to_scale = std(X(common_period),'omitmissing');

    % Normalize the reference
    mean_reference = mean(reference(common_period),'omitmissing');
    std_reference = std(reference(common_period),'omitmissing');
    
    scaled_index = (X - mean_to_scale) .* (std_reference ./ std_to_scale) + mean_reference;

end