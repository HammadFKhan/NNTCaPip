%% Parallel Batch processing 
function [directory, filename] = batchLoad()
directory = uigetdir(pwd,'Input Directory');
d = dir(directory);
for idx = 3:length(d)
   filename =  d(idx).name
   files = filename;
end

end

