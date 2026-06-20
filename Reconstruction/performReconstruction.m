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
        otherwise
            error('Invalid reconstruction method. Choose ''CPS'', ''PCR'', or ''PCA''.');
    end
end