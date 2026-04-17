% Run this code before anything else!
project_root = fileparts(mfilename('fullpath'));

addpath(genpath(project_root));
rmpath(genpath(fullfile(project_root, 'disciplines_old')));

fprintf("--- Paths Set Up ---\n");