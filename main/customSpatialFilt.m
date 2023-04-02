function PNR = customSpatialFilt(Y,params)
%%
use_sum = false;
if nargin < 4 || isempty(ROI_list)
    user_ROIs = 0;
else
    user_ROIs = 1;
    K = size(ROI_list,1);
end

dimY = ndims(Y) - 1;  % dimensionality of imaged data (2d or 3d)
sizY = size(Y);
T = sizY(end);        % # of timesteps
dx = sizY(1:dimY);    % # of voxels in each axis
d = prod(dx);         % total # of voxels  


if ~exist('params', 'var') params = []; end

if ~isfield(params, 'gSig') || isempty(params.gSig); 
    if dimY == 2; params.gSig = [5, 5]; else params.gSig = [5,5,5]; end
elseif length(params.gSig) == 1, params.gSig = params.gSig + zeros(1,dimY);
    if dimY == 3; params.gSig(3) = params.gSig(3)/2; end
end

if ~isfield(params, 'gSiz') || isempty(params.gSiz); 
    if ~iscell(params.gSig)
        params.gSiz = ceil(2*params.gSig + 1);
    else
        for j = 1:length(params.gSig)
            params.gSiz{j,1} = 2*params.gSig{j}+1; %cellfun(@times,params.gSig{j},num2cell(ones(size(params.gSig{j}))*2));
        end
    end
elseif length(params.gSiz) == 1, params.gSiz = params.gSiz + zeros(1,dimY);
    if dimY == 3; params.gSiz(3) = ceil(params.gSiz(3)/2); end
end

if isfield(params,'ssub'); 
    if ~iscell(params.gSig); params.gSig(1:2) = params.gSig(1:2)/params.ssub; params.gSiz(1:2) = ceil(params.gSiz(1:2)/params.ssub); 
    else
        for j = 1:length(params.gSig)
            params.gSig{j,1} = params.gSig{j}/params.ssub; %cellfun(@times,params.gSig{j},num2cell(ones(size(params.gSig{j}))/params.ssub));
            params.gSiz{j,1} = params.gSiz{j}/params.ssub; %cellfun(@times,params.gSiz{j},num2cell(ones(size(params.gSiz{j}))/params.ssub));
        end
    end
end

if ~isfield(params,'nb'), nb = 1; else nb = params.nb; end

if ~isfield(params, 'nIter'), nIter = 5; 
else nIter = params.nIter; end

if ~isfield(params, 'save_memory'), save_memory = 0;
else save_memory = params.save_memory; end
    
if ~isfield(params, 'chunkSiz'), chunkSiz = 100; else chunkSiz = params.chunkSiz; end
if ~isfield(params, 'windowSiz'), windowSiz = 32; else windowSiz = params.windowSiz; end
if ~isfield(params, 'med_app'), med_app = 1; else med_app = params.med_app; end

if ~isfield(params,'rem_prct') || isempty(params.rem_prct); params.rem_prct = 20; end
%%
psf = fspecial('gaussian', round(params.gSiz(1)), params.gSig(1));
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

% filter the data
HY = imfilter(Y, psf, 'replicate');
HY = reshape(HY, params.d1*params.d2, []);
% HY_med = median(HY, 2);
% HY_max = max(HY, [], 2)-HY_med;    % maximum projection
HY = bsxfun(@minus, HY, median(HY, 2));
HY_max = max(HY, [], 2);
Ysig = get_noise_fft(HY, params);
PNR = reshape(HY./Ysig, params.d1, params.d2,[]);
PNR0 = PNR;
PNR(PNR<params.min_pnr/5) = 0;