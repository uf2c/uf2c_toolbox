function varargout = uf2c_FuncInte(varargin)
% UF²C M-file for uf2c_FuncInte.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2019
%
% Copyright (c) 2019, Brunno Machado de Campos
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_FuncInte_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_FuncInte_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

warning('off','all')

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function uf2c_FuncInte_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_FuncInte_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function SSInp_Callback(hObject, eventdata, handles)
iROI = str2num(get(handles.SSInp,'String'));
sSize = 1+ iROI;
nVoxSeed = sSize^3;
VoxThre = round(nVoxSeed.*0.65);
set(handles.seedSizeTxt,'String', num2str(sSize))
set(handles.text33,'String', num2str(VoxThre))

function addFunc_Callback(hObject, eventdata, handles)
global fileFunc pathFunc numofVoluF nOFsubjects

set(handles.checkReg,'Enable','on')

try
    clear('regVets');
    delete([tmpDIR 'additional_Reg.mat']);
end

[fileFunc,pathFunc] = uigetfile({'*.nii','NIfTI files'},'Select all the functional images','MultiSelect','on');
if ~isequal(fileFunc,0)

    if ~iscell(fileFunc)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileFunc = {fileFunc};
    end
    fileFunc = sort(fileFunc); % SORT FILES IN THE ALPHABETIC ORDER
    nOFsubjects = size(fileFunc,2);

    for yt = 1:nOFsubjects
        tmpDir = dir([pathFunc fileFunc{yt}]);
        bitS(yt) = tmpDir.bytes;
    end

    sxs = unique(bitS);
    nDiftypes = size(sxs,2);
    if ~isequal(nDiftypes,1)
        warndlg({sprintf('There are %d distinct protocol types among the functional images that you added.',nDiftypes);...
          sprintf('Check the size of your files, it is the easiest way to identify the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
    end
    clear bitS

    set(handles.txtFunc,'String',sprintf('%d functional image(s) added',nOFsubjects))

    previewF = nifti([pathFunc fileFunc{1}]);
    matrixF = previewF.dat.dim(1:3);
    numofVoluF = previewF.dat.dim(4);
    VoxSizeF = previewF.hdr.pixdim(2:4);
    VoxS_Res = [num2str(VoxSizeF(1,1)),'x',num2str(VoxSizeF(1,2)),'x',num2str(VoxSizeF(1,3))];
    MatS_Res = [num2str(matrixF(1,1)),'x',num2str(matrixF(1,2)),'x',num2str(matrixF(1,3))];

    set(handles.numDyna,'String',numofVoluF)
    set(handles.SizeVoxFunc,'String',VoxS_Res)
    set(handles.SizeMatFunc,'String',MatS_Res)
    set(handles.txtregress,'String','0 added')
    set(handles.checkReg,'Value',0)
    set(handles.WSvalue,'String',numofVoluF)
    set(handles.pushbutton6,'Enable','on')
end

function addStru_Callback(hObject, eventdata, handles)
global fileStru pathStru

[fileStru,pathStru] = uigetfile({'*.nii','NIfTI files'},'Select all the Structural images','MultiSelect','on');
if ~isequal(fileStru,0)
    if ~iscell(fileStru)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileStru = {fileStru};
    end

    fileStru = sort(fileStru);

    for yt = 1:size(fileStru,2)
        tmpDir = dir([pathStru fileStru{yt}]);
        bitS(yt) = tmpDir.bytes;
    end
    sxs = unique(bitS);
    nDiftypes = size(sxs,2);
    if ~isequal(nDiftypes,1)
        warndlg({sprintf('There are %d distinct protocol types among the structural images that you added.',nDiftypes);...
          sprintf('Check the size of your files, it is the easiest way to identify the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
    end
    clear bitS

    set(handles.txtstru,'String',sprintf('%d Structural image(s) added',size(fileStru,2)))

    previewS = nifti([pathStru fileStru{1}]);
    matrixS = previewS.dat.dim;
    VoxSizeS = previewS.hdr.pixdim(2:4);
    VoxS_Res = [num2str(VoxSizeS(1,1)),'x',num2str(VoxSizeS(1,2)),'x',num2str(VoxSizeS(1,3))];
    MatS_Res = [num2str(matrixS(1,1)),'x',num2str(matrixS(1,2)),'x',num2str(matrixS(1,3))];

    set(handles.edit9,'String',VoxS_Res)
    set(handles.edit10,'String',MatS_Res)
    set(handles.pushbutton7,'Enable','on')
end

function Run_Callback(hObject, eventdata, handles)
global fileStru pathStru fileFunc pathFunc numofVoluF nrg tmpDIR

set(handles.status,'String','Running...')
drawnow
fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Process started\n');
fprintf('==========================================\r\n\r\n');

if isequal(get(handles.FazPrep,'Value'),0)
    if ~isequal(size(fileStru,2),size(fileFunc,2))
        warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
          'You need to add one functional image for each strutural and vice versa.';...
          'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
          'Attention: process aborted!');
        return
    end
end

SROI = str2num(get(handles.SSInp,'String')); % value to increase the ROI side
TR = str2num(get(handles.TRInp,'String'));    % repetition time (s)

try
    clear('regVets');
    delete([tmpDIR 'additional_Reg.mat']);
end

WS = str2num(get(handles.WSvalue,'String')); %Window Size (SIZE OF THE MOVING AVERAGE TO CALC DE CORRELATIONS)

ratioVol = numofVoluF/WS;

if ~isequal(ratioVol,round(ratioVol))
    warndlg('The "Window Size" and "Number of Dynamics" should to be multiples','Attention, runninng aborted!')
    return
end

multiWaitbar('Total Progress', 0, 'Color', 'g' );  %CHARGE WAIT BAR
if isequal(get(handles.FazPrep,'Value'),0)
    multiWaitbar('Filtering & Regression', 0, 'Color', 'y' ); %CHARGE WAIT BAR
end
multiWaitbar('Time series block', 0, 'Color', 'b');
multiWaitbar('Seed position', 0, 'Color', 'w' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

fideSJ = fopen([pathFunc,filesep,foldername,filesep,'1-Subjects.txt'],'w+'); 
for yTy = 1:size(fileFunc,2)
    fprintf(fideSJ,'%d: \t%s\r\n',yTy,fileFunc{1,yTy});
end
fclose(fideSJ);

imgRR = getframe(uf2c_FuncInte);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

Totmovtxt = fopen([pathFunc,filesep,foldername,filesep,'Total_Movement.txt'],'w+'); % CHARGE OUTPUT LOG FILE
fprintf(Totmovtxt,'Subject \t Maximum Displacemente (mm) \t Maximum Rotation (degree) \t Avg Framiwise Displacemente (FD) (mm) \t numb of FD Censored scans \t Avg DVARs (%%) \t numb of DVAR Censored scans \t Total Censored scans \r\n');


fideG = fopen([pathFunc,filesep,foldername,filesep,'Total_Log.txt'],'w+'); % CHARGE OUTPUT LOG FILE
fprintf(fideG,'                                          Functional Conectivity ROI Analysis\r\n\r\n');
fprintf(fideG,'Number of subjects included: %d \r\n\r\n', size(fileFunc,2));
fprintf(fideG,'Seed size: %dx%dx%d voxels³ \r\n\r\n', SROI+1, SROI+1, SROI+1);
fprintf(fideG, 'Subject_file_name                                       Number_of_Seeds         Summed_Corr_Value       Avg_Corr_Value      Summed_Right Corr_Value      Summed_Left Corr_Value      Avg_Right Corr_Value      Avg_Left Corr_Value\r\n');

if isequal(get(handles.FazPrep,'Value'),0) % just if you need preprocessing
%%%%%%%%%%%%%%%%%%%%%% FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isequal(get(handles.FazFilt,'Value'),1)

        Fs = 1/str2double(get(handles.TRInp,'String'));     % sampling frequency (Hz)

        FpassLOW = str2double(get(handles.LPF,'String'));   % passband frequency   (Hz) (default: 0.1)
        FstopLOW = str2double(get(handles.LPF,'String')) + (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.15)
        ApassLOW = 1;                                       % passband ripple      (dB) (default: 1)
        AstopLOW = str2double(get(handles.stba,'String'));  % stopband attenuation (dB) (default: 40)

        FstopHIGH = str2double(get(handles.HPF,'String')) - (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.005)
        FpassHIGH = str2double(get(handles.HPF,'String'));  % passband frequency   (Hz) (default: 0.008)
        AstopHIGH = str2double(get(handles.stba,'String')); % stopband attenuation (dB) (default: 40)
        ApassHIGH = 1;

        hLOW = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassLOW,FstopLOW,ApassLOW,AstopLOW,Fs);
        HdLOW = design(hLOW, 'equiripple');
        %fvtool(HdLOW)

        hHIGH  = fdesign.highpass('Fst,Fp,Ast,Ap',FstopHIGH,FpassHIGH,AstopHIGH,ApassHIGH,Fs);
        HdHIGH = design(hHIGH, 'equiripple');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

WinAvSavePos = struct;
WinAvSaveNeg = struct;

SPMdir2 = which('spm');
SPMdir2 = SPMdir2(1:end-5);

FunpreT = nifti([pathFunc,fileFunc{1}]);

if isequal(WS,FunpreT.dat.dim(4))
    dynStu= 0;
else
    dynStu= 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF SUBJECT LOOP %%%%%%%%%%%%%%%%%%%
for yy = 1: size(fileFunc,2)
   
    fprintf('\r\n')
    fprintf('UF²C ===================== %s\n',spm('time'));
    fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy}(1:14));
    fprintf('================================================\r\n\r\n');

    Funpre = nifti([pathFunc,fileFunc{yy}]);
    fvs = double(Funpre.hdr.pixdim(2:4)); %GET PIXEL SIZE (FUNCTIONAL)
    
    multiWaitbar('Total Progress', (yy/size(fileFunc,2))*0.95 );
    if isequal(get(handles.FazPrep,'Value'),0)
        multiWaitbar('Filtering & Regression', 'Reset' );
        dirname = Funpre.dat.fname;
        [ax,dirname,cx] = fileparts(dirname);
        
        if numel(dirname)>30
            fprintf('Attention!\n');
            fprintf('The filename is big. A short version will be used.\r\n')
            dirname = dirname(1:30);
        end
        nfidx = 1;
        while isequal(exist([pathFunc,dirname],'dir'),7)
            dirname = [dirname '_' num2str(nfidx)];
            nfidx = nfidx + 1;
        end
        
        mkdir(pathFunc,dirname)
        
        copyfile([pathStru,fileStru{yy}],[pathFunc,dirname,filesep,fileStru{yy}])
        file2 = fullfile(pathFunc,dirname,filesep,fileStru{yy});
        Strpre = nifti(file2);
        Svs = double(Strpre.hdr.pixdim(2:4)); % GET PIXEL SIZE (STRUCTURAL)
        copyfile([pathFunc,fileFunc{yy}],[pathFunc,dirname,filesep,fileFunc{yy}])
        
        file1 = cell(Funpre.dat.dim(4),1);
        for tt = 1:Funpre.dat.dim(4)
            file1{tt} = [pathFunc,dirname,filesep,fileFunc{yy},',',sprintf('%d',tt)];
        end
    else
        dirname = Funpre.dat.fname;
        [ax,dirname,cx] = fileparts(dirname);
        dirname = dirname(12:end);
        
        if numel(dirname)>30
            fprintf('Attention!\n');
            fprintf('The filename is big. A short version will be used.\r\n')
            dirname = dirname(1:30);
        end
        nfidx = 1;
        while isequal(exist([pathFunc,dirname],'dir'),7)
            dirname = [dirname '_' num2str(nfidx)];
            nfidx = nfidx + 1;
        end

        mkdir(pathFunc,dirname)
    end
    
    multiWaitbar('Time series block', 'Reset');
    multiWaitbar('Seed position', 'Reset');
    
    mkdir([pathFunc,dirname],'ROI_mask')
    mkdir([pathFunc,dirname],'Correlation_map')
    
    if isequal(get(handles.FazPrep,'Value'),0)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            STC = get(handles.STCo,'Value');
            BounB = 0;
            fileFuncX = fileFunc{yy};
            fileStruX = fileStru{yy};
            uf2c_Preproc(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX)
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d preprocessing. This is a SPM error. Check your data',yy), 'Proces Aborted')
            return
        end

        try
            close(fig);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot and save the motion series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mkdir([pathFunc,dirname,filesep,'Motion_Controls_Params'])
        Motfid = fopen([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'Motion_Quantifications.txt'],'w+');
        fprintf(Motfid,'Motion Quantifications\r\n\r\n');
        fprintf(Motfid,'Subject Name: \t%s\r\n',dirname);

        [Mot_fig,HDvV,HTvV] = uf2c_plot_motion([pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt'],'on');
        imgRR = getframe(Mot_fig);
        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot.png']);
        saveas(Mot_fig,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot'],'fig')

        fprintf(Motfid,'Maximum Displacemente (mm): %s \r\n',HDvV);
        fprintf(Motfid,'Maximum Rotation (degree): %s \r\n',HTvV);
        fprintf('Done! ============== %s\r\n',spm('time'));   
        
        fprintf(Totmovtxt,'%s \t %s \t',HDvV,HTvV);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BRAMILA's Framewise Displacemente (FD)
        % Code from: Brain and Mind Lab at Aalto University
        % Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also 
        % Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Quantifying Framewise Displacemente (FD)\r\n');

        if get(handles.FDcheck,'Value')
            cfg.motionparam = [pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt'];
            cfg.prepro_suite = 'spm';
            cfg.radius = 50;

            [FDts,rms] = bramila_framewiseDisplacement(cfg);

            save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement.mat'],'FDts')

            figuFD = figure('Visible','off');
            set(figuFD,'Name','Framewise Displacemente (FD) TS',...
                'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
                'Color',[1 0.94 0.86]);
            plot(FDts);
            hold on
            plot(ones(1,numel(FDts)).*str2num(get(handles.FDthres,'String')));
            ylabel('FD (mm)')
            xlabel('Time Points')
            title('Average Framewise Displacement TS','FontSize', 14);
            drawnow    

            imgRR = getframe(figuFD);
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS.png']);
            saveas(figuFD,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS'],'fig')
            close(figuFD)
            fprintf(Motfid,'BRAMILA''s average Framiwise Displacemente (FD): %.3f \r\n',mean(FDts));
            fprintf('Done! ============== %s\r\n',spm('time'));    
            fprintf(Totmovtxt,'%.5f \t',mean(FDts));
        else
            fprintf(Totmovtxt,'\t');
        end
        
        if get(handles.TM_FD,'Value')
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Creating FD Temporal Mask\r\n');

            FD_TM = FDts>str2num(get(handles.FDthres,'String'));
            fprintf('Done! ============== %s\r\n',spm('time'));
            fprintf(Totmovtxt,'%d \t',sum(FD_TM));
        else
            FD_TM = [];
            fprintf(Totmovtxt,'\t');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extracting Globals and Masking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            FI = 1;
            [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = uf2c_mask_globals(fvs,pathFunc,dirname,fileFunc{yy},fileStru{yy},0,FI);
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d globals and mask extraction. Check your data.',yy), 'Proces Aborted')
            return
        end
        try
            close(Mot_fig)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reshaping and Thresholding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            Dx = size(finalEPI,1);
            Dy = size(finalEPI,2);
            Dz = size(finalEPI,3);
            Dt = size(finalEPI,4);

            [finalEPI,histCut] = uf2c_imgthres(finalEPI,Dx,Dy,Dz,ScreSize,pathFunc,dirname);
            imgRR = getframe(histCut);
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Avg_Img_Histogram.tif']);
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d thresholding. Check your data.',yy), 'Proces Aborted')
            return
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BRAMILA DVARS - Computes Derivative VARiance
        % Code from: Brain and Mind Lab at Aalto University
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        extGMvar = get(handles.extGM,'Value');
        extWMvar = get(handles.extWM,'Value');
        extCSFvar = get(handles.extCSF,'Value');
        DVARthresV = str2num(get(handles.DVARthres,'String'));

        if extGMvar || extWMvar || extCSFvar
            [dvarsC1,dvarsC2,dvarsC3,MeanDvars] =  uf2c_bramila_dvars(pathFunc,dirname,finalEPI,extGMvar,extWMvar,extCSFvar,DVARthresV);
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) GM masked: %.3f \r\n',mean(dvarsC1));
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) WM masked: %.3f \r\n',mean(dvarsC2));
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) CSF masked: %.3f \r\n',mean(dvarsC3));
            fprintf(Totmovtxt,'%.4f \t',mean(MeanDvars));
        else
            fprintf(Totmovtxt,'\t');
        end

        if get(handles.TM_DVARs,'Value')
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Creating DVARs Temporal Mask\r\n');
            dvars_TM = zeros(size(finalEPI,4),1);
            dvars_TM = MeanDvars > str2num(get(handles.DVARthres,'String'));
            fprintf('Done! ============== %s\r\n',spm('time'));
            fprintf(Totmovtxt,'%d \t',sum(dvars_TM));
        else
            dvars_TM = [];
            fprintf(Totmovtxt,'\t');
        end

        fclose(Motfid);    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creating the regression Matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~get(handles.checkReg,'Value')
             nrg = 0;
        end

        if get(handles.checkMOV,'Value') || get(handles.checkOsc,'Value') || get(handles.checkReg,'Value')
            try
                [Frgval,Regstr,figu] = uf2c_regrprep(yy,get(handles.checkMOV,'Value'),pathFunc,dirname,fileFunc{yy},get(handles.checkOsc,'Value'),get(handles.checkReg,'Value'),ScreSize,MeanWM,MeanCSF,size(finalEPI,4),nrg);
                if ~isempty(figu)
                    if isequal(get(handles.TM_FD,'Value'),0) && isequal(get(handles.TM_DVARs,'Value'),0)
                        imgRR = getframe(figu);
                        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Regression Matrix.tif']);
                    end
                end
            catch
                set(handles.status,'String','Error...')
                warndlg(sprintf('An error occured during subject %d regressors creation. Check your data.',yy), 'Process  Aborted')
                return
            end
        else
            Regstr = '';
            Frgval = [];
        end            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adding the FD and DVAR regressors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TM_DVARsV = get(handles.TM_DVARs,'Value');
        TM_FDV = get(handles.TM_FD,'Value');

        if TM_DVARsV || TM_FDV
            [Frgval,figu,Censu] = uf2c_Add_FD_DVAR(figu,Frgval,pathFunc,dirname,TM_DVARsV,TM_FDV,FD_TM,dvars_TM,size(finalEPI,4));
            if size(Frgval,2)>1
                Regstr = [Regstr,'_tmpMask'];
            end
            nOfSensu = size(Censu,2);
            fprintf(Totmovtxt,'%d\r\n',nOfSensu);
        else
            fprintf(Totmovtxt,'\r\n',nOfSensu);
        end

        try
            close(figu)
        end
        try
            close(histCut)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Filtering and regressions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Performing regressions and/or filtering\r\n');

            BinFazFilt = get(handles.FazFilt,'Value');
            finalEPI = uf2c_filtregr(finalEPI,BinFazFilt,Dx,Dy,Dz,Dt,Frgval,HdLOW,HdHIGH,aFil,tmp1,tmp2,Regstr);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creating the FiltRegre Image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
                if  get(handles.FazFilt,'Value') && ~isequal(Regstr,'')
                    newEPI = Func;
                    newEPI.dat.fname = [pathFunc,dirname,filesep,'FI_FiltRegrSW_',fileFunc{yy}];
                    newEPI.descrip = 'UF²C-FI_Filtered_Regressed_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                else if get(handles.FazFilt,'Value') && isequal(Regstr,'')
                        newEPI = Func;
                        newEPI.dat.fname = [pathFunc,dirname,filesep,'FI_FiltSW_',fileFunc{yy}];
                        newEPI.descrip = 'UF²C-FI_Filtered_sw_EPI';
                        newEPI.dat.dtype = 'INT16-LE';
                        newEPI.dat(:,:,:,:) = finalEPI;
                        create(newEPI)
                        fprintf('Done! ============== %s\n\r',spm('time'));
                    else
                        newEPI = Func;
                        newEPI.dat.fname = [pathFunc,dirname,filesep,'FI_RegrSW_',fileFunc{yy}];
                        newEPI.descrip = 'UF²C-FI_Regressed_sw_EPI';
                        newEPI.dat.dtype = 'INT16-LE';
                        newEPI.dat(:,:,:,:) = finalEPI;
                        create(newEPI)
                        fprintf('Done! ============== %s\n\r',spm('time'));
                    end
                end
            else
                newEPI = Func;
                newEPI.dat.fname = [pathFunc,dirname,filesep,'FI_ThreSW_',fileFunc{yy}];
                newEPI.descrip = 'UF²C-FI_Masked and Thresholded_sw_EPI';
                newEPI.dat.dtype = 'INT16-LE';
                newEPI.dat(:,:,:,:) = finalEPI;
                create(newEPI)
                fprintf('Done! ============== %s\n\r',spm('time'));
            end
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d image creation. Check the avaiable disk space.',yy), 'Proces Aborted')
            return
        end
    else   % case preprocessing was not necessary
        Func = nifti([pathFunc,fileFunc{yy}]);
        finalEPI = Func.dat(:,:,:,:);
        Dx = size(finalEPI,1);
        Dy = size(finalEPI,2);
        Dz = size(finalEPI,3);
        Dt = size(finalEPI,4); 
    end

    fprintf('Each VOI will has a %d voxel side or %d voxels of volume \n', SROI+1, (SROI+1)^3)
    fprintf('Wait, this process can take several minutes \n')
    mt = zeros(Dx,Dy,Dz);
    nW=0;

    for i = 1:SROI+1:size(mt,1)-SROI
        for j = 1:SROI+1:size(mt,2)-SROI
            for k = 1:SROI+1:size(mt,3)-SROI
                try
                    mt = zeros(Dx,Dy,Dz);
                    mt(i,j,k) = 1;
                    mt(i,j:j+SROI,k) = 1;
                    mt(i:i+SROI,j,k) = 1;
                    mt(i:i+SROI,j:j+SROI,k) = 1;
                    mt(i,j,k:k+SROI) = 1;
                    mt(i,j:j+SROI,k:k+SROI) = 1;
                    mt(i:i+SROI,j,k:k+SROI) = 1;
                    mt(i:i+SROI,j:j+SROI,k:k+SROI) = 1;
                    VOIok = 0;
                catch
                     VOIok = 1;
                end
                if ~isequal(VOIok,1)
                    if nnz(mt.*finalEPI(:,:,:,1))>str2num(get(handles.text33,'String')) %threshold of minimum VOI size
                        nW = nW + 1;
                    end
                end
            end
        end
    end

    fprintf('There are %d possible positions to the cubic VOI. \n', nW)

    nro=0;
    finalEPI2cor = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';
    hg = 1;
    SEEDg = zeros(Dx,Dy,Dz);
    
    for MA = 1:WS:Dt
        multiWaitbar('Time series block', (MA+WS-1)/Dt);
        finalEPIx = finalEPI(:,:,:,MA:MA+(WS-1));
        finalEPI2corx = finalEPI2cor(MA:MA+(WS-1),:);
        Vn = 1;
        mapFt = zeros(Dx,Dy,Dz,nW);
        mapF = zeros(Dx,Dy,Dz,nW);
        for i = 1:SROI+1:size(mt,1)-SROI
            for j = 1:SROI+1:size(mt,2)-SROI
                for k = 1:SROI+1:size(mt,3)-SROI
                    multiWaitbar('Seed position',(Vn/(nW-1)));
                    try
                        mt = zeros(Dx,Dy,Dz);
                        mt(i,j,k) = 1;
                        mt(i,j:j+SROI,k) = 1;
                        mt(i:i+SROI,j,k) = 1;
                        mt(i:i+SROI,j:j+SROI,k) = 1;
                        mt(i,j,k:k+SROI) = 1;
                        mt(i,j:j+SROI,k:k+SROI) = 1;
                        mt(i:i+SROI,j,k:k+SROI) = 1;
                        mt(i:i+SROI,j:j+SROI,k:k+SROI) = 1;
                        VOIok = 0;
                    catch
                         VOIok = 1;
                    end
                    if ~isequal(VOIok,1)
                        if nnz(mt.*finalEPIx(:,:,:,1))> str2num(get(handles.text33,'String')) %threshold of minimum VOI size
                            if MA==1
                                SEEDg = SEEDg + (mt+Vn);
                            end
                            final4D = zeros(Dx,Dy,Dz,WS);
                            mean4Dnum = zeros(1,WS);
                            mean4Dden = zeros(1,WS);
                            for vx = 1:WS; % get the ROI voxels time series. This step exclude from the ROI eventuals with matter voxels.
                                final4D(:,:,:,vx) = mt.*finalEPIx(:,:,:,vx);
                                mean4Dnum(vx) = sum(sum(sum(final4D(:,:,:,vx))));
                                mean4Dden(vx) = nnz(final4D(:,:,:,vx));
                            end
                            mean4D = mean4Dnum./mean4Dden; %get the avegare ROI time series
                            map = corr(finalEPI2corx,mean4D');
                            map(isnan(map)) = 0;
                            mapFt(:,:,:,Vn) = reshape(map',Dx,Dy,Dz);
                            mapFt = squeeze(mapFt);
                            mapF(:,:,:,Vn) = mapFt(:,:,:,Vn) - (mt(:,:,:).*mapFt(:,:,:,Vn));
                            mapF2(:,:,:,Vn) = mapFt(:,:,:,Vn);
                            Vn = Vn+1;
                        else
                            nro = nro+1;
                        end
                    end
                end
            end
        end
        if MA == 1
            ROIstr = Funpre;
            ROIstr.dat.dim = [Dx Dy Dz];
            ROIstr.dat.fname = [pathFunc,dirname,filesep,'ROI_mask',filesep,'ROIsMask.nii'];
            ROIstr.dat.dtype = 'FLOAT32-LE';
            ROIstr.descrip = dirname;
            ROIstr.dat(:,:,:) = SEEDg;
            create(ROIstr)
        end
        asP = mapF>0;
        mapP = mapF.*asP;
        mapSOMA(:,:,:,hg) = sum(mapP,4);
        AVGmap(:,:,:,hg) = sum(abs(mapF),4)./(Vn-1);
        hg = hg+1;
    end
    save('mapF2.mat','mapF2','-v7.3');
    
    mapSOMA = squeeze(mapSOMA);
    AVGmap = squeeze(AVGmap);
    
    Mean4DVetPOS = zeros(1,(hg-1));
    Mean4DVetAVG = zeros(1,(hg-1));
    
    for ix = 1:(hg-1)
        baseP = mapSOMA(:,:,:,ix);
        NumCorWindowP = sum(sum(sum(baseP(baseP(:,:,:)>0))));
        DenCorWindowP = nnz(baseP(:,:,:)>0);
        corrWindowP = NumCorWindowP/DenCorWindowP;
        Mean4DVetPOS(ix) = corrWindowP;
        
        baseN = AVGmap(:,:,:,ix);
        NumCorWindowN = sum(sum(sum(baseN(baseN(:,:,:)>0))));
        DenCorWindowN = nnz(baseN(:,:,:)>0);
        corrWindowN = NumCorWindowN/DenCorWindowN;
        Mean4DVetAVG(ix) = corrWindowN;
    end
    
    meanPosR = mean(Mean4DVetPOS);
    meanAvgR = mean(Mean4DVetAVG);
    
    fig = figure;
    set(fig,'Name',dirname,'Position', [200 100 1350 800],...
                    'Color',[1 0.949019610881805 0.866666674613953]);
    subplot1 = subplot(1,2,1);
    plot(Mean4DVetPOS,'Parent',...
                    subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point' ,'FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Summed Positive Moving Window Average Cor value');
    
    subplot2 = subplot(1,2,2);
    plot(Mean4DVetAVG,'Parent',...
                    subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    
    title('Sum of positive values in each window (Interactivity)');
    
    Y1 = Mean4DVetPOS;
    Y2 = Mean4DVetAVG;

    box(subplot1,'on');
    hold(subplot1,'all');
    plot1 = plot(Y1,'Parent',subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.',...
        'LineWidth',1,...
        'DisplayName','data 1');
    xdata1 = get(plot1, 'xdata');
    ydata1 = get(plot1, 'ydata');
    xdata1 = xdata1(:);
    ydata1 = ydata1(:);
    nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
    if any(nanMask1)
        warning('GenerateMFile:IgnoringNaNs', ...
            'Data points with NaN coordinates will be ignored.');
        xdata1(nanMask1) = [];
        ydata1(nanMask1) = [];
    end
    axesLimits1 = xlim(subplot1);
    xplot1 = linspace(axesLimits1(1), axesLimits1(2));
    fitResults1 = polyfit(xdata1, ydata1, 2);
    yplot1 = polyval(fitResults1, xplot1);
    fitLine1 = plot(xplot1,yplot1,'DisplayName','quadratic',...
        'Parent',subplot1,'Tag','quadratic','Color',[0 0.75 0.75]);
    
    box(subplot2,'on');
    hold(subplot2,'all');
    plot2 = plot(Y2,'Parent',subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.',...
        'LineWidth',1,...
        'DisplayName','data 1');
    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Windows Averaged Corrrelation value (Interactivity)');
    xdata2 = get(plot2, 'xdata');
    ydata2 = get(plot2, 'ydata');
    xdata2 = xdata2(:);
    ydata2 = ydata2(:);
    nanMask2 = isnan(xdata2(:)) | isnan(ydata2(:));
    if any(nanMask1)
        warning('GenerateMFile:IgnoringNaNs', ...
            'Data points with NaN coordinates will be ignored.');
        xdata2(nanMask2) = [];
        ydata2(nanMask2) = [];
    end
    axesLimits2 = xlim(subplot2);
    xplot2 = linspace(axesLimits2(1), axesLimits2(2));
    fitResults2 = polyfit(xdata2, ydata2, 2);
    yplot2 = polyval(fitResults2, xplot2);
    fitLine2 = plot(xplot2,yplot2,'DisplayName','quadratic',...
    'Parent',subplot2,...
    'Tag','quadratic',...
    'Color',[0 0.75 0.75]);
    
    drawnow
    imgRR = getframe(gcf);
    imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Correlation_map',filesep,'Corr_Window_Graph.tif']);

    eval(sprintf('WinAvSavePos.suj_%d = Mean4DVetPOS;',yy));
    eval(sprintf('WinAvSaveAvg.suj_%d = Mean4DVetAVG;',yy));
        
    save([pathFunc,dirname,filesep,'Correlation_map',filesep,'Summed_Serie'],'Mean4DVetPOS')
    save([pathFunc,dirname,filesep,'Correlation_map',filesep,'Averaged_Serie'],'Mean4DVetAVG')

    WinCorr4DP = Func; 
    WinCorr4DA = Func; 
    
    WinCorr4DP.dat.dim = [Dx Dy Dz (hg-1)];
    WinCorr4DA.dat.dim = [Dx Dy Dz (hg-1)];

    WinCorr4DP.dat.fname = [pathFunc,dirname,filesep,'Correlation_map',filesep,dirname,'_Cr4D_Pos.nii'];
    WinCorr4DA.dat.fname = [pathFunc,dirname,filesep,'Correlation_map',filesep,dirname,'_Cr4D_Avg.nii'];
    
    WinCorr4DP.dat.dtype = 'FLOAT32-LE';
    WinCorr4DA.dat.dtype = 'FLOAT32-LE';
    
    WinCorr4DP.descrip = dirname;
    WinCorr4DA.descrip = dirname;

    WinCorr4DP.dat(:,:,:,1) = mapSOMA;
    WinCorr4DA.dat(:,:,:,1) = AVGmap;
    
    create(WinCorr4DP)
    create(WinCorr4DA)
    
    clear matlabbatch %Start to creates the smoothed correlation 4D map
    
    Corr4D_MultPOS = cell(1,(hg-1));
    
    for tx=1:(hg-1)
        Corr4D_MultPOS{tx} = [pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_Pos.nii' ',' num2str(tx)];
    end
    
    matlabbatch{1}.spm.spatial.smooth.data = Corr4D_MultPOS;                                       
    matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
    spm_jobman('run',matlabbatch)

    
    clear matlabbatch %Start to creates the smoothed correlation 4D map
    Corr4D_MultAVG = cell(1,(hg-1));
    
    for tx=1:(hg-1)
        Corr4D_MultAVG{tx} = [pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_Avg.nii' ',' num2str(tx)];
    end
    
    matlabbatch{1}.spm.spatial.smooth.data = Corr4D_MultAVG;                                       
    matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
    spm_jobman('run',matlabbatch)
    
    equa = 'i1';  %creates the equation to use in ImCalc
    if Dt>WS
        for eqn = 2:(Dt/WS)
            equa = sprintf('%s+i%d',equa,eqn);
        end
        equa = ['(' equa ')' sprintf('./%d',eqn)];
    
        Corr4D_MultP = transpose(Corr4D_MultPOS);
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = Corr4D_MultPOS;
        matlabbatch{1}.spm.util.imcalc.output = [dirname,'_MeanPOS.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc dirname filesep 'Correlation_map' filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = equa;
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        spm_jobman('run',matlabbatch)

        Corr4D_MultAVG = transpose(Corr4D_MultAVG);
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = Corr4D_MultAVG;
        matlabbatch{1}.spm.util.imcalc.output = [dirname,'_MeanAVG.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc dirname filesep 'Correlation_map' filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = equa;
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc dirname filesep 'Correlation_map' filesep dirname '_MeanPOS.nii']};                                       
        matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc dirname filesep 'Correlation_map' filesep dirname '_MeanAVG.nii']};                                       
        matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
        spm_jobman('run',matlabbatch)

        % positive values
            StruPF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_POS.nii']);
            matSF = StruPF.dat.dim;
            matlat = StruPF.dat(:,:,:,:);
                %Left values
                MaskR = zeros(matSF(1),matSF(2),matSF(3));
                MaskR(ceil(matSF(1)/2):matSF(1),:,:)=1;
                %Right values
                MaskL = zeros(matSF(1),matSF(2),matSF(3));
                MaskL(1:floor(matSF(1)/2),:,:)=1;
            leftR = zeros(1,matSF(4));
            rightR = zeros(1,matSF(4));

        % negative values
            StruNF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_AVG.nii']);   
            matSF2 = StruNF.dat.dim;
            matlat2 = StruNF.dat(:,:,:,:);
                %Left values
                MaskR2 = zeros(matSF2(1),matSF2(2),matSF2(3));
                MaskR2(ceil(matSF2(1)/2):matSF2(1),:,:)=1;
                %Right values
                MaskL2 = zeros(matSF2(1),matSF2(2),matSF2(3));
                MaskL2(1:floor(matSF2(1)/2),:,:)=1;
            leftR2 = zeros(1,matSF2(4));
            rightR2 = zeros(1,matSF2(4));

        for rt = 1:matSF(4)
            % positive values
                %Right values
                matlatL = matlat(:,:,:,rt).*MaskL;
                NumRLP = sum(sum(sum(matlatL)));
                DenRLP = nnz(matlatL(:,:,:));
                leftR(rt) = NumRLP/DenRLP;
                %Left values
                matlatR = matlat(:,:,:,rt).*MaskR;
                NumRRP = sum(sum(sum(matlatR)));
                DenRRP = nnz(matlatR(:,:,:));
                rightR(rt) = NumRRP/DenRRP;

            % negative values
                %Right values
                matlatL2 = matlat2(:,:,:,rt).*MaskL2;
                NumRLP2 = sum(sum(sum(matlatL2)));
                DenRLP2 = nnz(matlatL2(:,:,:));
                leftR2(rt) = NumRLP2/DenRLP2;
                %Left values
                matlatR2 = matlat2(:,:,:,rt).*MaskR2;
                NumRRP2 = sum(sum(sum(matlatR2)));
                DenRRP2 = nnz(matlatR2(:,:,:));
                rightR2(rt) = NumRRP2/DenRRP2;
        end
    
    else    %in case that there is no 4D volume (no window corr)
        
        % positive values
        StruPF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_POS.nii']);
        matSF = StruPF.dat.dim;
        matlat = StruPF.dat(:,:,:);
            %Left values
            MaskR = zeros(matSF(1),matSF(2),matSF(3));
            MaskR(ceil(matSF(1)/2):matSF(1),:,:)=1;
            %Right values
            MaskL = zeros(matSF(1),matSF(2),matSF(3));
            MaskL(1:floor(matSF(1)/2),:,:)=1;
    
        % negative values
            StruNF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_AVG.nii']);   
            matSF2 = StruNF.dat.dim;
            matlat2 = StruNF.dat(:,:,:,:);
            %Left values
            MaskR2 = zeros(matSF2(1),matSF2(2),matSF2(3));
            MaskR2(ceil(matSF2(1)/2):matSF2(1),:,:)=1;
            %Right values
            MaskL2 = zeros(matSF2(1),matSF2(2),matSF2(3));
            MaskL2(1:floor(matSF2(1)/2),:,:)=1;

            % positive values
                %Right values
                matlatL = matlat(:,:,:).*MaskL;
                NumRLP = sum(sum(sum(matlatL(matlatL(:,:,:)>0))));
                DenRLP = nnz(matlatL(:,:,:)>0);
                leftR = NumRLP/DenRLP;
                %Left values
                matlatR = matlat(:,:,:).*MaskR;
                NumRRP = sum(sum(sum(matlatR(matlatR(:,:,:)>0))));
                DenRRP = nnz(matlatR(:,:,:)>0);
                rightR = NumRRP/DenRRP;

            % negative values
                %Right values
                matlatL2 = matlat2(:,:,:).*MaskL2;
                NumRLP2 = sum(sum(sum(matlatL2(matlatL2(:,:,:)<0))));
                DenRLP2 = nnz(matlatL2(:,:,:)<0);
                leftR2 = NumRLP2/DenRLP2;
                %Left values
                matlatR2 = matlat2(:,:,:).*MaskR2;
                NumRRP2 = sum(sum(sum(matlatR2(matlatR2(:,:,:)<0))));
                DenRRP2 = nnz(matlatR2(:,:,:)<0);
                rightR2 = NumRRP2/DenRRP2;
    end
    
    if isequal(get(handles.OutputFil, 'Value'),1)
        delete([pathFunc dirname filesep fileFunc{yy}])
        delete([pathFunc dirname filesep fileStru{yy}])
        delete([pathFunc dirname filesep fileFunc{yy}(1:end-4) '.mat'])
        delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_inv_sn.mat'])
        delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_sn.mat'])
        delete([pathFunc dirname filesep 'mean' fileFunc{yy}])
        delete([pathFunc dirname filesep 'm' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc1' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc2' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc3' fileStru{yy}])
        delete([pathFunc dirname filesep 'w' fileFunc{yy}])
        delete([pathFunc dirname filesep 'w' fileFunc{yy}(1:end-4) '.mat'])
        delete([pathFunc dirname filesep 'sw' fileFunc{yy}(1:end-4) '.mat'])
        delete([pathFunc dirname filesep 'c1interp.nii'])
        delete([pathFunc dirname filesep 'sc1interp.nii'])
        delete([pathFunc dirname filesep 'c1interpF.nii'])
        delete([pathFunc dirname filesep 'c2interp.nii'])
        delete([pathFunc dirname filesep 'c3interp.nii'])
        delete([pathFunc dirname filesep 'Template.nii'])
        delete([pathFunc dirname filesep 'Template.mat'])
    end
    
    rightAverageP = mean(leftR);  % Yes, the name here should to be exchanged
    leftAverageP = mean(rightR);% Yes, the name here should to be exchanged
    
    rightAverageN = mean(leftR2);% Yes, the name here should to be exchanged
    leftAverageN = mean(rightR2);% Yes, the name here should to be exchanged

    fprintf(fideG,'%s:',dirname);% Write data in the output log file
    fprintf(fideG,'                           %d                       %.5f                            %.5f                           %.5f                            %.5f                           %.5f                            %.5f\r\n',nW,meanPosR,meanAvgR,rightAverageP,leftAverageP,rightAverageN,leftAverageN);% Write data in the output log file

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Subject %d, Done!\r\n',yy);

    
end  % End of the subject loop

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Generating overall results\r\n');


if size(fileFunc,2)>1
    for iy1 = 1:(hg-1)
        for ix3 = 1:size(fileFunc,2)
            eval(sprintf('STDPos_%d(%d) = WinAvSavePos.suj_%d(%d);',iy1,ix3,ix3,iy1));
            eval(sprintf('STDAvg_%d(%d) = WinAvSaveAvg.suj_%d(%d);',iy1,ix3,ix3,iy1));
        end
        StdFinalPos(iy1) = std(eval(sprintf('STDPos_%d',iy1)));
        StdFinalAvg(iy1) = std(eval(sprintf('STDAvg_%d',iy1)));
    end
    
    StdFinalPos = StdFinalPos./sqrt(size(fileFunc,2));
    StdFinalAvg = StdFinalAvg./sqrt(size(fileFunc,2));

end

graphfinalPos = WinAvSavePos.suj_1;
graphfinalAvg = WinAvSaveAvg.suj_1;
All_series_Pos(1,:) = WinAvSavePos.suj_1;
All_series_Avg(1,:) = WinAvSaveAvg.suj_1;

if size(fileFunc,2)>1
    for ix2 = 2:size(fileFunc,2)
        eval(sprintf('graphfinalPos = graphfinalPos + WinAvSavePos.suj_%d;',ix2));
        eval(sprintf('graphfinalAvg = graphfinalAvg + WinAvSaveAvg.suj_%d;',ix2));
        eval(sprintf('All_series_Pos(%d,:) = WinAvSavePos.suj_%d;',ix2,ix2));
        eval(sprintf('All_series_Avg(%d,:) = WinAvSaveAvg.suj_%d;',ix2,ix2));
    end
    
    graphfinalPos = graphfinalPos./size(fileFunc,2);
    graphfinalAvg = graphfinalAvg./size(fileFunc,2);
    
    save([pathFunc,filesep,foldername,filesep,'Pos_Average_Correlation_Serie'],'graphfinalPos')
    save([pathFunc,filesep,foldername,filesep,'Neg_Average_Correlation_Serie'],'graphfinalAvg')
    
    save([pathFunc,filesep,foldername,filesep,'All_Pos_Correlation_Series'],'All_series_Pos')
    save([pathFunc,filesep,foldername,filesep,'All_Avg_Correlation_Series'],'All_series_Avg')
end

try
    fig2 = figure;
    set(fig2,'Name','Total Moving Window Average','Position',[200 100 1350 800],...
                    'Color',[1 0.94 0.86]);
    
    subplot1 = subplot(1,2,1);
    plot(graphfinalPos,'Parent',...
                    subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point' ,'FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Summed interactivity');

    subplot2 = subplot(1,2,2);
    plot(graphfinalAvg,'Parent',...
                    subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Averaged interactivity');            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1 = graphfinalPos;
    Y2 = graphfinalAvg;

    box(subplot1,'on');
    hold(subplot1,'all');

    plot1 = errorbar(Y1,StdFinalPos/2,'Parent',subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',2,'Color',[0 0 1],'DisplayName','data 1');

    xdata1 = get(plot1, 'xdata');
    ydata1 = get(plot1, 'ydata');
    xdata1 = xdata1(:);
    ydata1 = ydata1(:);
    nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
    if any(nanMask1)
        warning('GenerateMFile:IgnoringNaNs', ...
            'Data points with NaN coordinates will be ignored.');
        xdata1(nanMask1) = [];
        ydata1(nanMask1) = [];
    end
    axesLimits1 = xlim(subplot1);
    xplot1 = linspace(axesLimits1(1), axesLimits1(2));
    fitResults1 = polyfit(xdata1, ydata1, 2);
    yplot1 = polyval(fitResults1, xplot1);
    fitLine1 = plot(xplot1,yplot1,'DisplayName','   quadratic',...
        'Parent',subplot1,...
        'Tag','quadratic',...
        'Color',[0 0.75 0.75]);
    box(subplot2,'on');
    hold(subplot2,'all');

    plot2 = errorbar(Y2,StdFinalAvg/2,'Parent',subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',2,'Color',[0 0 1],'DisplayName','data 1');

    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Averaged Moving Window Interactivity value');
    xdata2 = get(plot2, 'xdata');
    ydata2 = get(plot2, 'ydata');
    xdata2 = xdata2(:);
    ydata2 = ydata2(:);
    nanMask2 = isnan(xdata2(:)) | isnan(ydata2(:));
    if any(nanMask1)
        warning('GenerateMFile:IgnoringNaNs', ...
            'Data points with NaN coordinates will be ignored.');
        xdata2(nanMask2) = [];
        ydata2(nanMask2) = [];
    end
    axesLimits2 = xlim(subplot2);
    xplot2 = linspace(axesLimits2(1), axesLimits2(2));
    fitResults2 = polyfit(xdata2, ydata2, 2);
    yplot2 = polyval(fitResults2, xplot2);
    fitLine2 = plot(xplot2,yplot2,'DisplayName','   quadratic',...
        'Parent',subplot2,...
        'Tag','quadratic',...
        'Color',[0 0.75 0.75]);

    saveas(fig2,[pathFunc,foldername,filesep,'Totals_Windows_Graph'],'jpg')
catch
    
end
set(handles.status,'String','Done!')

multiWaitbar('CloseAll'); %close the progressbar window
fclose('all');
fprintf('\r\n')
fprintf('All Done! ========== %s\n',spm('time'));
fprintf('==========================================\r\n');


function addreg_Callback(hObject, eventdata, handles)
global regVets
try
    clear(regVets)
end
uf2c_Addreg

function checkReg_Callback(hObject, eventdata, handles)
if isequal(get(handles.checkReg,'Value'),1)
    set(handles.addreg,'Enable','on')
    set(handles.refreshb,'Enable','on')
else
    set(handles.addreg,'Enable','off')
    set(handles.refreshb,'Enable','off')
    set(handles.txtregress,'String','0 added')
end

function refreshb_Callback(hObject, eventdata, handles)
global nrg 
try
    if nrg < 0
        nrg = 0;
    end
    nrT = num2str(nrg);
    set(handles.txtregress,'String',[nrT ' added'])
end

if nrg == 0
     set(handles.checkReg,'Value',0)
     set(handles.addreg,'Enable','off')
     set(handles.refreshb,'Enable','off')
end

function FazFilt_Callback(hObject, eventdata, handles)
if isequal(get(handles.FazFilt,'Value'),1)
    set(handles.LPF,'Enable','on')
    set(handles.HPF,'Enable','on')
    set(handles.stba,'Enable','on')
else
    set(handles.LPF,'Enable','off')
    set(handles.HPF,'Enable','off')
    set(handles.stba,'Enable','off')
end

function FazPrep_Callback(hObject, eventdata, handles)
if get(handles.FazPrep,'Value')
    set(handles.addStru,'Enable','off')
    set(handles.addFunc,'String','Add FiltRegrSW Files')
    set(handles.FazFilt,'Enable','off')
    set(handles.LPF,'Enable','off')
    set(handles.HPF,'Enable','off')
    set(handles.stba,'Enable','off')
    set(handles.checkMOV,'Enable','off')
    set(handles.checkOsc,'Enable','off')
    set(handles.checkReg,'Enable','off')
    set(handles.text14,'Enable','off')
    set(handles.text15,'Enable','off')
    set(handles.text21,'Enable','off')
    set(handles.HPFtxt,'Enable','off')
    set(handles.LPFtxt,'Enable','off')
    set(handles.DVAR_TS,'Enable','off')
else
    set(handles.addStru,'Enable','on')
    set(handles.addFunc,'String','Add Functional Files')
    set(handles.FazFilt,'Enable','on')
    set(handles.LPF,'Enable','on')
    set(handles.HPF,'Enable','on')
    set(handles.stba,'Enable','on')
    set(handles.checkMOV,'Enable','on')
    set(handles.checkOsc,'Enable','on')
    set(handles.checkReg,'Enable','on')
    set(handles.text14,'Enable','on')
    set(handles.text15,'Enable','on')
    set(handles.text21,'Enable','on')
    set(handles.HPFtxt,'Enable','on')
    set(handles.LPFtxt,'Enable','on')
    set(handles.DVAR_TS,'Enable','off')
end

function SPMbb_Callback(hObject, eventdata, handles)
if get(handles.SPMbb,'Value')
    set(handles.SSInp,'String','3')
    set(handles.seedSizeTxt,'String','4')
    set(handles.text33,'String','42')
else
    set(handles.SSInp,'String','5')
    set(handles.seedSizeTxt,'String','6')
    set(handles.text33,'String','140')
end

function pushbutton6_Callback(hObject, eventdata, handles)
global fileFunc

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Func.txt']);
end

subL = transpose(fileFunc);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Func.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileFunc{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Func.txt'])

function pushbutton7_Callback(hObject, eventdata, handles)
global fileStru

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Stru.txt']);
end

subL = transpose(fileStru);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Stru.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileStru{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Stru.txt'])

function WSvalue_Callback(hObject, eventdata, handles)
function OutputFil_Callback(hObject, eventdata, handles)
function checkMOV_Callback(hObject, eventdata, handles)
function stba_Callback(hObject, eventdata, handles)
function checkOsc_Callback(hObject, eventdata, handles)
function HPF_Callback(hObject, eventdata, handles)
function TRInp_Callback(hObject, eventdata, handles)
function LPF_Callback(hObject, eventdata, handles)
function SizeMatFunc_Callback(hObject, eventdata, handles)
function edit10_Callback(hObject, eventdata, handles)
function edit9_Callback(hObject, eventdata, handles)
function SizeVoxFunc_Callback(hObject, eventdata, handles)
function numDyna_Callback(hObject, eventdata, handles)

function WSvalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TRInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function stba_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function numDyna_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SizeVoxFunc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SizeMatFunc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SSInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Saveori4D_Callback(hObject, eventdata, handles)
