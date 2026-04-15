% Run this code before anything else!
project_root = fileparts(mfilename('fullpath'));

addpath(genpath(project_root));
rmpath(genpath(fullfile(project_root, 'disciplines_old')));

fprintf('Path configured. Project root:\n  %s\n\n', project_root);
fprintf('Folders on path:\n');
fprintf('  constraints/\n');
fprintf('  data/\n');
fprintf('  disciplines/\n');
fprintf('  tests/\n');
fprintf('  simplified_problem/\n');