function h = nanimagesc(varargin)
% a wrapper for imagesc, with some formatting going on for nans

% plotting data. Removing and scaling axes (this is for image plotting)
h = imagesc(varargin{:});

[~,idim] = max(cellfun(@getdims, varargin));
img_data = varargin{idim};

% setting alpha values
if ismatrix(img_data)
  set(h, 'AlphaData', ~isnan(img_data))
elseif ndims(img_data) == 3
  set(h, 'AlphaData', ~isnan(img_data(:, :, 1)))
end

if nargout < 1
  clear h
end