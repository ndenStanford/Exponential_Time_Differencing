% clear output in debug_plot
% Specify the folder where the files live.
myFolder = 'C:\Users\nden\Desktop\old stuffs\Stanford\Research\ENERGY223\Project\project4_Nutchapol_Dendumrongsup_IMPES\debug_plots';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  % fprintf(1, 'Now deleting %s\n', fullFileName);
  delete(fullFileName);
end