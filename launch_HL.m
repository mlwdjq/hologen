
% Add this root directory to path
[cDirThis, cName, cExt] = fileparts(mfilename('fullpath'));
addpath(genpath(cDirThis));
addpath('../mpm');

% mic library
mpm addpath


