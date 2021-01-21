%% this script is used to generate Quadriwave zone plates
if ~exist('hl')||~ishandle(hl.hFigure)
    launch_HL;
end

%% parameter setting (SHARP)
hl.uieLambda.set(13.5);
hl.uieFocalLength.set(0.5);
hl.uieNA.set(0.0825);
hl.uieOffsetAngle.set(6);
hl.uieObscuration.set(0);
hl.uieMinFeature.set(1);
hl.uieOffset.set('[0, 0.05255]');
hl.uieShift.set('[0, 0]');
hl.uipHologram.setSelectedIndex(uint8(1));% QWLSI
hl.cb(hl.uipHologram);

%% file path, select gds file here
savePath = fullfile(hl.cAppPath,'..','..','Data','gds','Quadriwave zone plates');
if ~isdir(savePath)
    mkdir(savePath);
end
hl.uieFilePath.set(savePath);

%% loop
T_um = [135,90,45,30,22.5,18,13.5,9,6.75,4.5,2.7,1.8,1.35];
for i = 1:length(T_um)
    hl.uieT.set(T_um(i));
    hl.uieFileName.set(['onAxis_13.5nm_0.0825_',num2str(T_um(i)),'um']);
    
    %% simulate pattern
    hl.cb(hl.uibGenPattern);
    drawnow;
    
    %% simulate pattern
    hl.cb(hl.uibGenGDS);
end