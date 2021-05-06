
function movie_data = FastTiff
% [filename, pathname] = uigetfile({'*.tiff;*.tif'}, 'Pick a image video file');
% if isequal(filename,0) || isequal(pathname,0)
%    disp('User pressed cancel')
% else
%    disp(['User selected ', fullfile(pathname, filename)])
% end
H = waitbar(0,'Compressing Tiff Files');
warning('off','all') % Suppress all the tiff warnings
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
movie_data = zeros(I,J,K);
movie_data(:,:,1)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    movie_data(:,:,n) = tstack.read();
end
warning('on','all')
delete(H)
