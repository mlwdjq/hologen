
% Add this root directory to path
[cDirThis, cName, cExt] = fileparts(mfilename('fullpath'));
addpath(genpath(cDirThis));
addpath('../mpm');

% mic library
mpm addpath


%% generating holo lens for 1d LSI (visible light)
% genHL_LSI_onAxis_Visible();

%% generating holo lens for 2d LSI (visible light)
% genHL_QWLSI_onAxis_Visible();

%% generating holo lens for 2d LSI (EUV light)
 genHL_QWLSI_onAxis_EUV();
