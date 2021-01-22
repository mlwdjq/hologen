%% test script
if ~exist('hl')||~ishandle(hl.hFigure)
    launch_HL;
end

%% parameter setting
hl.uieLambda.set(13.5);
hl.uieT.set(1.35);
hl.uieFocalLength.set(0.5);
hl.uieNA.set(0.0225);
hl.uieOffsetAngle.set(6);
hl.uieObscuration.set(0);
hl.uieMinFeature.set(1);
hl.uieFileName.set('test');
hl.uieOffset.set('[0, 0]');
hl.uieShift.set('[0, 0]');
hl.uipHologram.setSelectedIndex(uint8(1));
hl.cb(hl.uipHologram);

%% file path, select gds file here
savePath = fullfile(hl.cAppPath,'..','..','Data','gds');
if ~isdir(savePath)
    mkdir(savePath);
end
hl.uieFilePath.set(savePath);


%% simulate pattern
hl.cb(hl.uibGenPattern);
drawnow;

%% simulate pattern
hl.cb(hl.uibGenGDS);
