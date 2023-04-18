function [Img_p,fileInfo] = CaImAnReadIn(fileName)

% ScanImageTiffReader for Matlab
%
% This is a Matlab package for extracting data from
% <https://en.wikipedia.org/wiki/Tagged_Image_File_Format Tiff> and
% <http://bigtiff.org/ BigTiff> files recorded using
% <http://scanimage.org ScanImage>.  It is a very fast tiff reader and provides
% access to ScanImage-specific metadata.  It should read most tiff files, but as
% of now we don't support compressed or tiled data.  It is also available as a
% <https://vidriotech.gitlab.io/scanimagetiffreader-julia/ Julia>,
% <https://vidriotech.gitlab.io/scanimagetiffreader-python/ Python>,  or 
% <https://vidriotech.gitlab.io/scanimage-tiff-reader C library>.  There's also a
% <https://vidriotech.gitlab.io/scanimage-tiff-reader command line interface>.
%
% Both <http://scanimage.org ScanImage> and this reader are products of 
% <https://vidriotechnologies.com/ Vidrio Technologies>.  If you have
% questions or need support feel free to <https://gitlab.com/vidriotech/scanimagetiffreader-matlab/issues submit an issue>
% <https://vidriotechnologies.com/contact-support contact us>. 
%
%%% Downloads
% Packages includes mex files built against 64-bit Matlab 2016b (v9.1)
% and targeting Windows, Linux, and OS X. 
%
% <html>
% <table>
%   <tr><td><b>Target</b></td><td></td><td><b>Version</b></td></td>
%   <tr><td>Windows x64</td><td><a href="https://gitlab.com/api/v4/projects/8668010/jobs/artifacts/master/download?job=build_windows">Download</a></td><td>v1.3</td></tr>
%   <tr><td>OS X/Linux</td><td>Please contact support@mbfbioscience.com with your request!</td><td>v1.3</td></tr>
% </table>
% </html>
%
%%% Installation
%
% # Download the build for your operating system.
% # Unzip it.
% # Copy the |+ScanImageTiffReader| folder to a location on your Matlab path.
%
%% Examples
% Import the reader class.  A constructed reader represents an open file.
import ScanImageTiffReader.ScanImageTiffReader;
%% 
% *Read a stack*
%
% Note that Matlab's ordering of dimensions means the image might be 
% transposed from what you expect.  We leave the transpose up to you.
fprintf('Reading Scanimage file...');
reader=ScanImageTiffReader(fileName);
vol=reader.data();
vol = permute(vol,[2 1 3]);
fprintf('done\n');
fprintf('Applying Kalman filter...')
Img_p= kalman_stack_filter(vol);
fprintf('done\n')

%%
if plotOn
    figure,
    set(gcf, 'Units', 'pixels')
    set(gcf, 'Position',[100,100, 1200, 520]/1.2)
    for i = 1:5:size(vol,3)
        subplot(1,2,1),imagesc(vol(:,:,i)),colormap(gray), title(['Frame: ' num2str(i)])
        subplot(1,2,2),imagesc(Img_p(:,:,i)),colormap(gray), title(['Frame: ' num2str(i)]),caxis([min(Img_p(:,:,i),[],'all') max(Img_p(:,:,i),[],'all')])
        pause(0.033);
    end
end
%% 
% *Get some metadata!*
meta=reader.metadata();
desc=reader.descriptions();
fileInfo.desc = desc;
fileInfo.meta = meta;

