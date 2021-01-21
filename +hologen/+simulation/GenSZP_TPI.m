%% this script is used to generate shearing zone plates
if ~exist('hl')||~ishandle(hl.hFigure)
    launch_HL;
end

%% parameter setting (TPI)
hl.uieLambda.set(13.5);
hl.uieFocalLength.set(3);
hl.uieNA.set(0.0825);
hl.uieOffsetAngle.set(6);
hl.uieObscuration.set(0.025);
hl.uieMinFeature.set(1);
hl.uieOffset.set('[0, 0]');
hl.uieShift.set('[0, 0]');
hl.uipHologram.setSelectedIndex(uint8(2));% LSI
hl.cb(hl.uipHologram);

%% file path, select gds file here
savePath = fullfile(hl.cAppPath,'..','..','Data','gds','Shearing zone plates');
if ~isdir(savePath)
    mkdir(savePath);
end
hl.uieFilePath.set(savePath);

%% loop
T_um = [810,540,270,180,135,108,81,54,40.5,27,16.2,10.8,8.1];
for i = 1:length(T_um)
    hl.uieT.set(T_um(i));
    hl.uieFileName.set(['onAxis_13.5nm_0.0825_',num2str(T_um(i)),'um_1D']);
    
    %% simulate pattern
    hl.cb(hl.uibGenPattern);
    drawnow;
    
    %% simulate pattern
    hl.cb(hl.uibGenGDS);
end
