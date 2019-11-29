function varargout = uf2c_JustPreprocHypo(varargin)
% UF²C M-file for uf2c_JustPreprocHypo.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2017
%
% Copyright (c) 2017, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @uf2c_JustPreprocHypo_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_JustPreprocHypo_OutputFcn, ...
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

function uf2c_JustPreprocHypo_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_JustPreprocHypo_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function AddFunc_Callback(hObject, eventdata, handles)
global fileFunc pathFunc numofVoluF nOFsubjects matrixF VoxSizeF

set(handles.checkReg,'Enable','on')

try
    clear('regVets');
end
try
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
    set(handles.checFunc,'Enable','on')
end

function AddStru_Callback(hObject, eventdata, handles)
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
    drawnow
    previewS = nifti([pathStru fileStru{1}]);
    matrixS = previewS.dat.dim;
    VoxSizeS = previewS.hdr.pixdim(2:4);

    VoxS_Res = [num2str(VoxSizeS(1,1)),'x',num2str(VoxSizeS(1,2)),'x',num2str(VoxSizeS(1,3))];
    MatS_Res = [num2str(matrixS(1,1)),'x',num2str(matrixS(1,2)),'x',num2str(matrixS(1,3))];
    set(handles.edit9,'String',num2str(VoxS_Res))
    set(handles.edit10,'String',num2str(MatS_Res))
    set(handles.Run,'Enable','on')
    set(handles.checStru,'Enable','on')
end

function Run_Callback(hObject, eventdata, handles)
global fileStru pathStru fileFunc pathFunc numofVoluF nrg tmpDIR

if ~isequal(size(fileStru,2),size(fileFunc,2))
    warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
      'You need to add one functional image for each strutural and vice versa.';...
      'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
      'Attention: process aborted!');
    return
end

set(handles.status,'String','Running....')
drawnow

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Process  started\n');
fprintf('==========================================\r\n');


TR = str2num(get(handles.TRInp,'String'));    % repetition time (s)

multiWaitbar('Total Progress', 0, 'Color', 'g' );  %CHARGE WAIT BAR
multiWaitbar('Filtering & Regression', 0, 'Color', 'y' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('1-Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

fideSJ = fopen([pathFunc,filesep,foldername,filesep,'1-Subjects.txt'],'w+'); 
for yTy = 1:size(fileFunc,2)
    fprintf(fideSJ,'%d: \t%s\r\n',yTy,fileFunc{1,yTy});
end
fclose(fideSJ);

imgRR = getframe(JustPreprocHypo);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

Totmovtxt = fopen([pathFunc,filesep,foldername,filesep,'Total_Movement.txt'],'w+'); % CHARGE OUTPUT LOG FILE
fprintf(Totmovtxt,'Subject \t Maximum Displacemente (mm) \t Maximum Rotation (degree) \t Avg Framiwise Displacemente (FD) (mm) \t numb of FD Censored scans \t Avg DVARs (%%) \t numb of DVAR Censored scans \t Total Censored scans \r\n');

if get(handles.FazFilt,'Value') % just if you want band-pass filtering
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Designing filters\n');
    %%%%%%%%%%%%%%%%%%%%%% FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fs = 1/TR;     % sampling frequency (Hz)
    aFil = 1;
    
    % LowPass design
    FpassLOW = str2double(get(handles.LPF,'String'));   % passband frequency   (Hz) (default: 0.1)
    if ~isnan(FpassLOW)
        FstopLOW = str2double(get(handles.LPF,'String')) + (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.15)
        ApassLOW = 1;                                       % passband ripple      (dB) (default: 1)
        AstopLOW = str2double(get(handles.stba,'String'));  % stopband attenuation (dB) (default: 40)
        
        hLOW = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassLOW,FstopLOW,ApassLOW,AstopLOW,Fs);
        HdLOW = design(hLOW, 'equiripple');
        tmp1 = 3*(max(length(HdLOW.Numerator),length(aFil))-1);
    else
        HdLOW = [];
        tmp1 = [];
    end
    
    % HighPass design
    FpassHIGH = str2double(get(handles.HPF,'String'));  % passband frequency   (Hz) (default: 0.008)
    if ~isnan(FpassHIGH)
        FstopHIGH = str2double(get(handles.HPF,'String')) - (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.005)
        if FstopHIGH<0
            FstopHIGH = str2double(get(handles.HPF,'String'))-(str2double(get(handles.HPF,'String')).*0.9);
        end
        AstopHIGH = str2double(get(handles.stba,'String')); % stopband attenuation (dB) (default: 40)
        ApassHIGH = 1;

        hHIGH  = fdesign.highpass('Fst,Fp,Ast,Ap',FstopHIGH,FpassHIGH,AstopHIGH,ApassHIGH,Fs);
        HdHIGH = design(hHIGH, 'equiripple');
        tmp2 = 3*(max(length(HdHIGH.Numerator),length(aFil))-1);
    else
        HdHIGH = [];
        tmp2 = [];
    end
    
    fprintf('Done! ============== %s\r\n',spm('time'));
else
    FpassLOW = [];  % passband frequency   (Hz) (default: 0.1)
    FstopLOW = []; % stopband frequency   (Hz) (default: 0.15)
    ApassLOW = [];                                     % passband ripple      (dB) (default: 1)
    AstopLOW = [];  % stopband attenuation (dB) (default: 40)
    FpassHIGH = [];  % passband frequency   (Hz) (default: 0.008)
    FstopHIGH = []; % stopband frequency   (Hz) (default: 0.005)
    AstopHIGH = []; % stopband attenuation (dB) (default: 40)
    ApassHIGH = [];
    hLOW = [];
    HdLOW = [];
    hHIGH  = [];
    HdHIGH = [];
    aFil = [];
    tmp1 = [];
    tmp2 = [];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

WinAvSave = struct;
ExclSub = 0;

SPMdir2 = which('spm');
SPMdir2 = SPMdir2(1:end-5);

for yy = 1: size(fileFunc,2) % START OF SUBJECT LOOP
    
    multiWaitbar('Total Progress', (yy/size(fileFunc,2))*0.95 );
    multiWaitbar('Filtering & Regression', 'Reset' );

    fprintf('\r\n')
    fprintf('UF²C ===================== %s\n',spm('time'));
    if size(fileFunc{yy},2)<15
        fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy});
    else
        fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy}(1:14));
    end
    fprintf('================================================\r\n');
    
    try
        close(fig)
    end
    
    Funpre = nifti([pathFunc,fileFunc{yy}]);
    fvs = double(Funpre.hdr.pixdim(2:4)); %GET PIXEL SIZE (FUNCTIONAL)
    
    fvs(3) = fvs(2);
    
    
    dirname = Funpre.dat.fname;
    [ax,dirname,cx] = fileparts(dirname);
    
    if numel(dirname)>30
        fprintf('Attention!\n');
        fprintf('The filename is long. A shorter version will be used.\r\n')
        dirname = dirname(1:30);
    end
    
    nfidx = 1;
    while isequal(exist([pathFunc,dirname],'dir'),7)
        dirname = [dirname '_' num2str(nfidx)];
        nfidx = nfidx + 1;
    end
    
    fprintf(Totmovtxt,'%s\t',dirname);
    
    mkdir(pathFunc,dirname)

    copyfile([pathStru,fileStru{yy}],[pathFunc,dirname,filesep,fileStru{yy}])
    file2 = fullfile(pathFunc,dirname,filesep,fileStru{yy});
    Strpre = nifti(file2);
    Svs = double(Strpre.hdr.pixdim(2:4)); % GET PIXEL SIZE (STRUCTURAL)
    copyfile([pathFunc,fileFunc{yy}],[pathFunc,dirname,filesep,fileFunc{yy}])
    
    file1 = cell(Funpre.dat.dim(4),1);
    for tt = 1:Funpre.dat.dim(4)
        file1{tt,1} = [pathFunc,dirname,filesep,fileFunc{yy},',',sprintf('%d',tt)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
         STC = get(handles.STCo,'Value');
         BounB = get(handles.SPMbb,'Value');
         fileFuncX = fileFunc{yy};
         fileStruX = fileStru{yy};
         Preproc_uf2c_Hypota(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX)
    catch
         set(handles.status,'String','Error...')
         warndlg(sprintf('An error occured during subject %d preprocessing. This is a SPM error. Check your data',yy), 'Process  Aborted')
         return
    end
    try
        close(fig);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot and save the motion series
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkdir([pathFunc,dirname,filesep,'Motion_Controls_Params'])
    Motfid = fopen([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'1-Motion_Quantifications.txt'],'w+');
    fprintf(Motfid,'Motion Quantifications\r\n\r\n');
    
    
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
        fprintf(Totmovtxt,'N/A\t');
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
        fprintf(Totmovtxt,'N/A\t');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extracting Globals and Masking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        FI = 0;
        [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = mask_globals_uf2c(fvs,pathFunc,dirname,fileFunc{yy},fileStru{yy},get(handles.SPMbb,'Value'),FI);
    catch
        set(handles.status,'String','Error...')
        warndlg(sprintf('An error occured during subject %d globals and mask extraction. Check your data.',yy), 'Process  Aborted')
        return
    end
    try
        close(Mot_fig)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshapeing and Thresholding
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        Dx = size(finalEPI,1);
        Dy = size(finalEPI,2);
        Dz = size(finalEPI,3);
        Dt = size(finalEPI,4);
        
        [finalEPI,histCut] = imgthres_uf2c(finalEPI,Dx,Dy,Dz,ScreSize,pathFunc,dirname);
        imgRR = getframe(histCut);
        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Avg_Img_Histogram.tif']);
    catch
        set(handles.status,'String','Error...')
        warndlg(sprintf('An error occured during subject %d thresholding. Check your data.',yy), 'Process  Aborted')
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
        [dvarsC1,dvarsC2,dvarsC3,MeanDvars] = bramila_dvars_uf2c(pathFunc,dirname,finalEPI,extGMvar,extWMvar,extCSFvar,DVARthresV);
        fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) GM masked: %.3f \r\n',mean(dvarsC1));
        fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) WM masked: %.3f \r\n',mean(dvarsC2));
        fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) CSF masked: %.3f \r\n',mean(dvarsC3));
        fprintf(Totmovtxt,'%.4f \t',mean(MeanDvars));
    else
        MeanDvars = [];
        fprintf(Totmovtxt,'N/A\t');
    end
    
    if get(handles.TM_DVARs,'Value')
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Creating DVARs Temporal Mask\r\n');
        dvars_TM = zeros(size(finalEPI,4),1);
        dvars_TM = MeanDvars>str2num(get(handles.DVARthres,'String'));
        fprintf('Done! ============== %s\r\n',spm('time'));
        fprintf(Totmovtxt,'%d \t',sum(dvars_TM));
    else
        dvars_TM = [];
        fprintf(Totmovtxt,'N/A\t');
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
            [Frgval,Regstr,figu,thicksL] = regrprep_uf2c(yy,get(handles.checkMOV,'Value'),pathFunc,dirname,fileFunc{yy},get(handles.checkOsc,'Value'),get(handles.checkReg,'Value'),ScreSize,MeanWM,MeanCSF,size(finalEPI,4),nrg);
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
        thicksL = {};
    end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding the FD and DVAR regressors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TM_DVARsV = get(handles.TM_DVARs,'Value');
    TM_FDV = get(handles.TM_FD,'Value');
    DVAR_tsV = get(handles.DVAR_TS,'Value');
    
    if TM_DVARsV || TM_FDV || DVAR_tsV
        [Frgval,figu,Censu] = Add_FD_DVAR_uf2c(figu,thicksL,Frgval,pathFunc,dirname,TM_DVARsV,TM_FDV,DVAR_tsV,MeanDvars,FD_TM,dvars_TM,size(finalEPI,4));
        if size(Frgval,2)>1
            Regstr = [Regstr,'_tmpMask'];
        end
        if TM_DVARsV || TM_FDV
            nOfSensu = size(Censu,2);
            fprintf(Totmovtxt,'%d\r\n',nOfSensu);
        end
    else
        fprintf(Totmovtxt,'N/A\r\n');
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
    if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'') || TM_DVARsV || TM_FDV
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Performing regressions and/or filtering\r\n');

        BinFazFilt = get(handles.FazFilt,'Value');
        finalEPI = filtregr_uf2cHypo(finalEPI,BinFazFilt,Dx,Dy,Dz,Dt,Frgval,HdLOW,HdHIGH,aFil,tmp1,tmp2,Regstr);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creating the FiltRegre Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
            if  get(handles.FazFilt,'Value') && ~isequal(Regstr,'')
                newEPI = Func;
                newEPI.dat.fname = [pathFunc,dirname,filesep,'FiltRegrSW_',fileFunc{yy}];
                newEPI.descrip = 'UF²C-Filtered_Regressed_sw_EPI';
                newEPI.dat.dtype = 'INT16-LE';
                newEPI.dat(:,:,:,:) = finalEPI;
                create(newEPI)
                fprintf('Done! ============== %s\n\r',spm('time'));
            else if get(handles.FazFilt,'Value') && isequal(Regstr,'')
                    newEPI = Func;
                    newEPI.dat.fname = [pathFunc,dirname,filesep,'FiltSW_',fileFunc{yy}];
                    newEPI.descrip = 'UF²C-Filtered_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                else
                    newEPI = Func;
                    newEPI.dat.fname = [pathFunc,dirname,filesep,'RegrSW_',fileFunc{yy}];
                    newEPI.descrip = 'UF²C-Regressed_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                end
            end
        else
            newEPI = Func;
            newEPI.dat.fname = [pathFunc,dirname,filesep,'ThreSW_',fileFunc{yy}];
            newEPI.descrip = 'UF²C-Masked and Thresholded_sw_EPI';
            newEPI.dat.dtype = 'INT16-LE';
            newEPI.dat(:,:,:,:) = finalEPI;
            create(newEPI)
            fprintf('Done! ============== %s\n\r',spm('time'));
        end
    catch
        set(handles.status,'String','Error...')
        warndlg(sprintf('An error occured during subject %d image creation. Check the avaiable disk space.',yy), 'Process   Aborted')
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deleting residual images and files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isequal(get(handles.OutputFil, 'Value'),1)
        delete([pathFunc dirname filesep fileStru{yy}])
        delete([pathFunc dirname filesep fileFunc{yy}])
        delete([pathFunc dirname filesep fileFunc{yy}(1:end-4) '.mat'])
        delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg8.mat'])
        delete([pathFunc dirname filesep 'mean' fileFunc{yy}])
        delete([pathFunc dirname filesep 'y_' fileStru{yy}])
        delete([pathFunc dirname filesep 'm' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc1' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc2' fileStru{yy}])
        delete([pathFunc dirname filesep 'wc3' fileStru{yy}])
        
        if get(handles.STCo,'Value')
            delete([pathFunc dirname filesep 'wa' fileFunc{yy}])
            delete([pathFunc dirname filesep 'a' fileFunc{yy}])
        end
        
        delete([pathFunc dirname filesep 'sw' fileFunc{yy}])
        delete([pathFunc dirname filesep 'c1interp.nii'])
        delete([pathFunc dirname filesep 'sc1interp.nii'])
        delete([pathFunc dirname filesep 'c1interpF.nii'])
        delete([pathFunc dirname filesep 'c2interp.nii'])
        delete([pathFunc dirname filesep 'c3interp.nii'])
        delete([pathFunc dirname filesep 'Template.nii'])
        delete([pathFunc dirname filesep 'Template.mat'])
    end
    
    clear newEPI finalEPIresh Frgval rgval
    fprintf('Subject %d Done! ====== %s\n\r',yy,spm('time'));
    
end

multiWaitbar('CloseAll'); %close the progressbar window
set(handles.status,'String','Done!!!')

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('All Done!\n');
fprintf('==========================================\r\n');

    
function checkReg_Callback(hObject, eventdata, handles)
if isequal(get(handles.checkReg,'Value'),1)
    set(handles.addreg,'Enable','on')
    set(handles.refreshb,'Enable','on')
else
    set(handles.addreg,'Enable','off')
    set(handles.refreshb,'Enable','off')
    set(handles.txtregress,'String','0 added')
end

function addreg_Callback(hObject, eventdata, handles)
global regVets

try
    clear(regVets)
end
Addreg

function refreshb_Callback(hObject, eventdata, handles)
global nrg 
try
    if nrg < 0
        nrg=0;
    end
    nrT = num2str(nrg);
    set(handles.txtregress,'String',[nrT ' added'])
end

if nrg == 0
     set(handles.checkReg,'Value',0)
     set(handles.addreg,'Enable','off')
     set(handles.refreshb,'Enable','off')
end

function checFunc_Callback(hObject, eventdata, handles)
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

function checStru_Callback(hObject, eventdata, handles)
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

function stba_Callback(hObject, eventdata, handles)
function FazFilt_Callback(hObject, eventdata, handles)
function LPF_Callback(hObject, eventdata, handles)
function SPMbb_Callback(hObject, eventdata, handles)
function OutputFil_Callback(hObject, eventdata, handles)
function checkMOV_Callback(hObject, eventdata, handles)
function checkOsc_Callback(hObject, eventdata, handles)
    
function HPF_Callback(hObject, eventdata, handles)
    
function TRInp_Callback(hObject, eventdata, handles)
function SavReg_Callback(hObject, eventdata, handles)
function STCo_Callback(hObject, eventdata, handles)
function DVARthres_Callback(hObject, eventdata, handles)
function FDthres_Callback(hObject, eventdata, handles)

function HPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function stba_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TRInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FDthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DVARthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extGM_Callback(hObject, eventdata, handles)
if get(handles.extGM,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extWM,'Value') || get(handles.extCSF,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end
    
function extWM_Callback(hObject, eventdata, handles)
if get(handles.extWM,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extCSF,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end

function extCSF_Callback(hObject, eventdata, handles)
if get(handles.extCSF,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
else
    if get(handles.extGM,'Value') || get(handles.extWM,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end

function FDcheck_Callback(hObject, eventdata, handles)
if get(handles.FDcheck,'Value')
    set(handles.TM_FD,'Value',0)
    set(handles.TM_FD,'Enable','on')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.TM_FD,'Value',0)
    set(handles.TM_FD,'Enable','off')
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
end

function TM_FD_Callback(hObject, eventdata, handles)
if get(handles.TM_FD,'Value')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
end

function TM_DVARs_Callback(hObject, eventdata, handles)
if get(handles.TM_DVARs,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
else
    set(handles.text27,'Enable','off')
    set(handles.DVARthres,'Enable','off')
end

function DVAR_TS_Callback(hObject, eventdata, handles)
