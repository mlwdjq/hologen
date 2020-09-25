classdef HL_generator < mic.Base
    
    
    properties (Constant)
        dWidth  = 870;
        dHeight =  600;
        
        % Axes tab IDS
        U8SIM          = 1
    end
    
    properties
        
        cAppPath = fileparts(mfilename('fullpath'))
        
        % Graphical elements
        hFigure     % Main figure (not overwritable)
        
        % valuables
        dGDSData
        dGDSDataHead
        dGDSDataEnd
        dPolygon
        dMask
        dFarfield
        dPupil
        dPupilAmp
        dPhase
        dRMS
        dZrn
        dUx
        dUy
        x_um
        y_um
        yp_um
        xp_um
        
        % axis tab
        uitgAxesDisplay     % displays axes:simulation
        haFieldAmp
        

        % parameters
        hpPara
        uieLambda
        uieT
        uieFocalLength
        uieNA
        uieOffsetAngle
        uieObscuration
        uieMinFeature
        uieFileName
        uieOffset
        uieShift
        uipHologram
        
        % control
        hpControl
        uieFilePath
        uibSetPath
        uilFileList
        uibGenGDS
        uibOpenGDS
        uibGenPattern
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        function this = HL_generator()
            this.init()
        end
        
        
        
        function init(this)
            
            % axis tab
            this.uitgAxesDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Simulation'});
           
            % parameters
            this.uieLambda       = mic.ui.common.Edit('cLabel', 'Wavelength(nm)', 'cType', 'd');
            this.uieT     = mic.ui.common.Edit('cLabel', 'T(nm)', 'cType', 'd');
            this.uieFocalLength = mic.ui.common.Edit('cLabel', 'Focal length(mm)', 'cType', 'd');
            this.uieNA  = mic.ui.common.Edit('cLabel', 'NA', 'cType', 'd');
            this.uieOffsetAngle  = mic.ui.common.Edit('cLabel', 'Offset angle(deg)', 'cType', 'd');
            this.uieObscuration      = mic.ui.common.Edit('cLabel', 'Obs. radius(mm)', 'cType', 'd');
            this.uieMinFeature      = mic.ui.common.Edit('cLabel', 'Min feature(nm)', 'cType', 'd');
            this.uieFileName      = mic.ui.common.Edit('cLabel', 'File name', 'cType', 'c');
            this.uieOffset         = mic.ui.common.Edit('cLabel', 'Center offset(mm)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieShift         = mic.ui.common.Edit('cLabel', 'Center shift(mm)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipHologram    = mic.ui.common.Popup('cLabel', 'Hologram', 'ceOptions',...
                {'QWLSI','LSI'}, 'lShowLabel', true);
            this.uieLambda.set(13.5);
            this.uieT.set(54);
            this.uieFocalLength.set(3);
            this.uieNA.set(0.0875);
            this.uieOffsetAngle.set(6);
            this.uieObscuration.set(0);
            this.uieMinFeature.set(1);
            this.uieFileName.set('New');
            this.uieOffset.set('[0, 0]');
            this.uieShift.set('[0, 0]');
            this.uipHologram.setSelectedIndex(uint8(1));
            
            % control
            this.uieFilePath    = mic.ui.common.Edit('cLabel', 'GDS file path', 'cType', 'c','fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSetPath    = mic.ui.common.Button('cText', 'Set path', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uilFileList         = mic.ui.common.List('cLabel', 'GDS file list', ...
                'lShowDelete', false, 'lShowMove', false, 'lShowRefresh', false);
            this.uibGenGDS   = mic.ui.common.Button('cText', 'Gen GDS', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibGenPattern   = mic.ui.common.Button('cText', 'Gen pattern', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibOpenGDS   = mic.ui.common.Button('cText', 'Open GDS', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieFilePath.set('');
            

        end
        
        % Callback handler
        function cb(this, src,evt)
            switch src
                case {this.uieOffset, this.uieShift}
                    this.validateCouplesEditBox(src, '[-1, 1]');
                    
                case this.uieFilePath
                    filename = this.uieFilePath.get();
                    if ~isempty(filename)
                        if ~strcmp(filename,'Please load GDS file first!')
                            this.dataLoading(filename);
                        end
                    else
                        this.uieFilePath.set('Please load GDS file first!');
                    end
                    
                case this.uibSetPath
                    cDataDir = fullfile(this.cAppPath, '..','..','data', '*.gds');
                    [d, p] = uigetfile(cDataDir);
                    if isequal(d,0)||isequal(p,0)
                        return;
                    end
                    filename = [p d];
                    this.uieFilePath.set(filename);
                    this.dataLoading(filename);
                    
                case this.uibGenGDS
                    this.genGDS();
                    
                case this.uibGenPattern
                    this.genPattern();
                    
                case this.uibOpenGDS
                    this.openGDS();
                    

            end
        end
        
        function dataLoading(this, filename)
            % load data
            tic,
            outputFile=fopen(filename,'rb');
            if outputFile==-1
                fprintf('Opening file failed!\n');
                return;
            end
            offset=0;
            while 1
                fseek(outputFile, offset, 'bof');
                temp=fread(outputFile, [4,1],'uint8');
                if temp(3)==8&&temp(4)==0
                    break;
                else
                    offset=temp(1)*16+temp(2)+offset;
                end
            end
            fseek(outputFile, offset, 'bof');
            gds=fread(outputFile, 'int32','b');
            gds(end-1:end)=[];
            s1=find(gds==264192);% head of the polygon
            s2=find(gds==266496);% end of the polygon
            % remove error datas
            ds1 = diff(s1);
            ds2 = diff(s2);
            for i = 1:length(ds1)
                try
                    while ds1(i)>ds2(i)
                        ds2(i) = ds2(i)+ds2(i+1);
                        ds2(i+1) =[];
                        s2(i+1)  =[];
                    end
                    while ds1(i)<ds2(i)
                        ds1(i) = ds1(i)+ds1(i+1);
                        ds1(i+1) =[];
                        s1(i+1)  =[];
                    end
                catch
                    break;
                end
            end
            num=length(s1);
            fclose(outputFile);
            if num~=length(s2)
                fprintf('Data format is not correct!\n');
                return;
            end
            this.dGDSData = gds;
            this.dGDSDataHead = s1;
            this.dGDSDataEnd = s2;
            fprintf('Reading GDS file took %ds\n',round(toc));
        end
        
        % propagate function
        function propagate(this)
            % --  Exposure tool information
            s2 = this.dGDSDataEnd;
            s1 = this.dGDSDataHead;
            gds = this.dGDSData;
            if isempty(gds)
                fprintf('Please reload datafile!\n');
                return;
            end
            tic,
            num=length(s1);
            wavl = this.uieLambda.get();
            nrd = this.uieT.get();
            offsetAngle = this.uieOffsetAngle.get()*pi/180;
            azimuth = this.uieObscuration.get()*pi/180;
            propdis_nm=this.uieFocalLength.get()*1e6;
            defocus_mm = this.uieNA.get();
            propMethod = this.uipHologram.getSelectedIndex();
            if propMethod ==1
                sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
                    (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2)+...% point source illumination
                    2*pi/wavl*sqrt(x.^2+y.^2+...
                    (propdis_nm-x*tan(offsetAngle)*cos(azimuth)-y*tan(offsetAngle)*sin(azimuth)).^2)); % compensate phase;
            else
                sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
                    (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2));% point source illumination
            end
            rangeX = eval(this.uieOffset.get());
            rangeY = eval(this.uieShift.get());
            xLeft= rangeX(1);
            xRight= rangeX(2);
            yLeft= rangeY(1);
            yRight= rangeY(2);
            
            xc=linspace(sin(atan(xLeft/propdis_nm*1000)),sin(atan(xRight/propdis_nm*1000)),2*nrd-1)/wavl;
            yc=linspace(sin(atan(yLeft/propdis_nm*1000)),sin(atan(yRight/propdis_nm*1000)),2*nrd-1)/wavl;
            x_um=propdis_nm*tan(asin(xc*wavl))/1000;
            y_um=propdis_nm*tan(asin(yc*wavl))/1000;
            [x_nm,y_nm] = meshgrid(x_um*1000,y_um*1000);
            % [x_nm,y_nm] = meshgrid(linspace(xLeft,xRight,2*nrd-1)*1000,linspace(yLeft,yRight,2*nrd-1)*1000);
            [nyux,nyuy]=meshgrid(xc,yc); % frequency coordinates
            this.dUx = nyux;
            this.dUy = nyuy;
            this.x_um = x_um;
            this.y_um = y_um;
            dpx_um = wavl*propdis_nm/(xRight-xLeft)*1e-6;
            dpy_um = wavl*propdis_nm/(yRight-yLeft)*1e-6;
            this.xp_um = dpx_um*[-nrd+1:nrd-1];
            this.yp_um = dpy_um*[-nrd+1:nrd-1];
            %% create polygon structs
            polyg(num,1)=struct('xy',[],'tx',[],'phase',[]);
            backg.int=0;
            backg.phase=0;
            for j=1:num
                if abs(cos(azimuth))==1
                    xs=gds((s1(j)+5):2:(s2(j)-3))/10*cos(offsetAngle); % from A to nm (tilt only works for x axis)
                    ys=gds((s1(j)+6):2:(s2(j)-3))/10;
                else
                    xs=gds((s1(j)+5):2:(s2(j)-3))/10; % from A to nm (tilt only works for x axis)
                    ys=gds((s1(j)+6):2:(s2(j)-3))/10*cos(offsetAngle);
                end
                polyg(j).xy =[xs';ys'];
                polyg(j).tx = 1;
                x0=mean(xs);
                y0=mean(ys);
                polyg(j).phase = sph(x0,y0);
                %     figure(2),plot3(x0,y0,sph(x0,y0),'.'),hold on;
            end
            
            % Forces all polygons to be defined in a CW orientation
            polyg = GDS.utils.cc_or_ccw(polyg) ;
            
            % -- Pupil intensity calculation
            % Computes the fourier transform of polygons in the plane specified by
            if propMethod ==1
                this.dFarfield = GDS.utils.polyProp(polyg, backg, nrd,nyux,nyuy) ;
            else
                this.dFarfield = GDS.utils.sphereProp(polyg,x_nm,y_nm,propdis_nm,wavl,offsetAngle,azimuth);
            end
            this.dPolygon = polyg;
            fprintf('Propagation took %ds\n',round(toc));
            % Make Field tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8SIM);
            this.replot(this.U8SIM);
        end
        
        
        
        
        
        
        % validates whether a char edit box evaluates to a Nx2 matrix,
        % colors accordingly.  Empty value is changed to []
        function [lOut, vals] = validateCouplesEditBox(~, src, cDefaultVal)
            lOut = true;
            vals = [];
            if isempty(src.get())
                src.styleDefault();
                src.set(cDefaultVal);
                return
            end
            try
                vals = eval(src.get());
                [~, sc] = size(vals);
                if (sc == 2 || sc == 0)
                    src.styleDefault();
                    lOut = true;
                else
                    src.styleBad();
                    lOut = false;
                end
            catch
                % can't read this edit box
                src.styleBad();
                lOut = false;
            end
        end
        
        % Main redraw function. Pass tab indices to refresh axes
        function replot(this, dTabIdx)
            
            switch dTabIdx
                
                case this.U8SIM
                    diff_amp = this.dFarfield;
                    norm_flag = 1 ; % 1: Normalize the intensity, 0: Unnormalize the intensity
                    if norm_flag>0
                        diff_amp = diff_amp./max(abs(diff_amp(:))) ;
                    else
                        diff_amp = diff_amp./length(diff_amp).^2 ;
                    end
                    fieldAmp=abs(diff_amp);
                    fieldPha=atan2(imag(diff_amp),real(diff_amp));
                    
                    imagesc(this.haFieldAmp,this.x_um,this.y_um,fieldAmp);axis(this.haFieldAmp,'xy'); colorbar(this.haFieldAmp);
                    xlabel(this.haFieldAmp,'x/um'),ylabel(this.haFieldAmp,'y/um');
                   
               
            end
            
        end
        
        
        
        
        function build(this, hFigure, dOffsetX, dOffsetY)
            if nargin <3
                dOffsetX = 0;
                dOffsetY = 0;
            elseif nargin == 3
                dOffsetY = dOffsetX;
                dOffsetX =  hFigure;            
            end
            
            % build the main window
            if nargin == 2||nargin == 4
                this.hFigure = hFigure;
            else
                this.hFigure = figure(...
                    'name', 'GDS propagation GUI v1.200701',...
                    'Units', 'pixels',...
                    'Position', [5 - dOffsetX, 5 - dOffsetY,  this.dWidth, this.dHeight],...
                    'handlevisibility','off',... %out of reach gcf
                    'numberTitle','off',...
                    'Toolbar','none',...
                    'Menubar','none');
            end
            
            
            % Build all containers first:
            drawnow
            
            % Axes
            dTgPx = 20;
            dTgPy = 20;
            this.uitgAxesDisplay.build(this.hFigure, dTgPx, dTgPy, 550, 550);
            
            % Axes:Far field
            uitField = this.uitgAxesDisplay.getTabByName('Simulation');

            
            this.haFieldAmp = axes('Parent', uitField, ...
                'Units', 'pixels', ...
                'Position', [50, 60, 430, 360], ...
                'XTick', [], 'YTick', []);

          

            
            % parameters
            this.hpPara = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Parameters',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 310 260 270] ...
                );
            this.uieLambda.build (this.hpPara,20,20,100,20);
            this.uieT.build (this.hpPara,140,20,100,20);
            this.uieFocalLength.build (this.hpPara,20,60,100,20);
            this.uieNA.build (this.hpPara,140,60,100,20);
            this.uieOffsetAngle.build (this.hpPara,20,100,100,20);
            this.uieObscuration.build (this.hpPara,140,100,100,20);
            this.uieMinFeature.build (this.hpPara,20,140,100,20);
            this.uieFileName.build (this.hpPara,140,140,100,20);
            this.uieOffset.build (this.hpPara,20,180,100,20);
            this.uieShift.build (this.hpPara,140,180,100,20);
            this.uipHologram.build (this.hpPara,20,220,100,20);
            
            % control
            this.hpControl = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Control',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 30 260 270] ...
                );
            this.uieFilePath.build (this.hpControl,20,20,220,20);
            this.uibSetPath.build (this.hpControl,140,60,100,20);
            this.uilFileList.build (this.hpControl,20,83,220,105);
            this.uibGenGDS.build (this.hpControl,140,240,100,20);
            this.uibGenPattern.build (this.hpControl,20,240,100,20);
            this.uibOpenGDS.build (this.hpControl,140,210,100,20);
            
            drawnow;
        end
        
    end
    
end
