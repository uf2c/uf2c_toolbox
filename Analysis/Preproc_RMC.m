function varargout = Preproc_RMC(varargin)
% UF²C M-file for Preproc_RMC.fig
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
                   'gui_OpeningFcn', @Preproc_RMC_OpeningFcn, ...
                   'gui_OutputFcn',  @Preproc_RMC_OutputFcn, ...
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

function Preproc_RMC_OpeningFcn(hObject, eventdata, handles, varargin)
global PsOri ps
if license('test','Distrib_Computing_Toolbox')
    ps = parallel.Settings;
    PsOri = ps.Pool.AutoCreate;
end

handles.output = hObject;
guidata(hObject, handles);

function varargout = Preproc_RMC_OutputFcn(hObject, eventdata, handles) 
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
    set(handles.LVR,'String',floor(numofVoluF/10))
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
global fileStru pathStru fileFunc pathFunc numofVoluF nrg tmpDIR PsOri ps

UDV = uf2c_defaults('Preproc_RMC');

EPCo = get(handles.EPCo,'Value');

if EPCo
    vxx = ver;
    delete(gcp('nocreate'))
    if license('test','Distrib_Computing_Toolbox')
        numcores = feature('numcores');
        try
            parpool(numcores-UDV.NC2R);
        end
    else
        set(handles.EPCo,'Value',0);
    end
else
    if license('test','Distrib_Computing_Toolbox')
        ps.Pool.AutoCreate = false;
    end
end

multiWaitbar('CloseAll'); %close the progressbar window

if ~isequal(size(fileStru,2),size(fileFunc,2))
    warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
      'You need to add one functional image for each strutural and vice versa.';...
      'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
      'Attention: process aborted!');
   fclose('all');
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
multiWaitbar('Total Preprocessing Progress', 0, 'Color', 'y' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('1-Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

fideSJ = fopen([pathFunc,filesep,foldername,filesep,'1-Subjects.txt'],'w+'); 
for yTy = 1:size(fileFunc,2)
    fprintf(fideSJ,'%d: \t%s\r\n',yTy,fileFunc{1,yTy});
end
fclose(fideSJ);

imgRR = getframe(Preproc_RMC);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

options = struct('format','pdf','outputDir',[pathFunc,filesep,foldername,filesep]);
publish('uf2c_defaults.m',options);

Totmovtxt = fopen([pathFunc,filesep,foldername,filesep,'Total_Movement.txt'],'w+'); % CHARGE OUTPUT LOG FILE
fprintf(Totmovtxt,'Subject \t Maximum Displacemente (mm) \t Maximum Rotation (degree) \t Avg Framiwise Displacemente (FD) (mm) \t numb of FD suprathres scans \t Avg DVARs (%%) \t numb of DVAR suprathres scans \t Total suprathres scans \r\n');

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
multiWaitbar('Total Progress', 0.25);

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Imaging Processing and Quantifications\r\n');

for yy = 1: size(fileFunc,2) % START OF SUBJECT LOOP
    
    multiWaitbar('Total Preprocessing Progress', (yy/size(fileFunc,2)).*0.95);
    
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
    dirname = Funpre.dat.fname;
    [ax,dirname,cx] = fileparts(dirname);
    
    if numel(dirname)>UDV.MSDN
        fprintf('Attention!\n');
        fprintf('The filename is long. A shorter version will be used.\r\n')
        dirname = dirname(1:UDV.MSDN);
    end
    
    nfidx = 1;
    while isequal(exist([pathFunc,dirname],'dir'),7)
        dirname = [dirname '_' num2str(nfidx)];
        nfidx = nfidx + 1;
    end
    
    mkdir([pathFunc,dirname])
    fprintf(Totmovtxt,'%s \t',dirname);

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
         Preproc_uf2c(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX)
    catch
         set(handles.status,'String','Error...')
         warndlg(sprintf('An error occured during subject %d preprocessing. This is a SPM error. Check your data',yy), 'Process  Aborted')
        fclose('all');
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
    
        [Mot_fig,HDvV,HTvV] = uf2c_plot_motion([pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt'],UDV.MFVstr);
    
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
        cfg.radius = UDV.FDBR;

        [FDts,~] = bramila_framewiseDisplacement(cfg);

        save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement.mat'],'FDts')

        figuFD = figure('Visible','off');
        set(figuFD,'Name','Framewise Displacemente (FD) TS',...
            'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
            'Color',[1 0.94 0.86]);
        
        BinVetFDMEAN = double(FDts>str2num(get(handles.FDthres,'String')));
        BinVetFDMEAN = BinVetFDMEAN.*FDts;

        plot(FDts);
        hold on
        plot(ones(1,numel(FDts)).*str2num(get(handles.FDthres,'String')));
        area(BinVetFDMEAN,'FaceColor','r')
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
        
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Creating FD Temporal Mask\r\n');

        FD_TM = FDts>str2num(get(handles.FDthres,'String'));
        fprintf('Done! ============== %s\r\n',spm('time'));
        fprintf(Totmovtxt,'%d \t',sum(FD_TM));
        fprintf(Motfid,'Number of FD suprathreshold volumes: %d \r\n',sum(FD_TM));
        DVAR_FD_BIN(1) = 1;
    else
        fprintf(Totmovtxt,'N/A\t');
        FD_TM = [];
        DVAR_FD_BIN(1) = 0;
        fprintf(Totmovtxt,'N/A\t');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extracting Globals and Masking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        FI = UDV.FI;
        [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = mask_globals_uf2c(fvs,pathFunc,dirname,fileFunc{yy},fileStru{yy},get(handles.SPMbb,'Value'),FI,UDV.WMI,UDV.NativeS);
    catch
        set(handles.status,'String','Error...')
        warndlg(sprintf('An error occured during subject %d globals and mask extraction. Check your data.',yy), 'Process  Aborted')
       fclose('all');
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
       fclose('all');
        return
    end
    
    %%%%%% Smothing
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Smoothing\r\n');
    
    newEPI = Func;
    newEPI.dat.fname = [pathFunc,dirname,filesep,'tmpW_',fileFunc{yy}];
    newEPI.descrip = 'UF²C-Masked and Thresholded_sw_EPI';
    newEPI.dat.dtype = 'INT16-LE';
    newEPI.dat(:,:,:,:) = finalEPI;
    create(newEPI)

    clear file1
    for tt = 1:newEPI.dat.dim(4)
        file1{tt,1} = [pathFunc,dirname,filesep,'tmpW_',fileFunc{yy},',',sprintf('%d',tt)];
    end
    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = file1;
    if BounB
        matlabbatch{1}.spm.spatial.smooth.fwhm = [UDV.SmoFac*fvs(1) UDV.SmoFac*fvs(2) UDV.SmoFac*fvs(3)];
    else
        matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    end
    
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',matlabbatch)
    
    preFiReStru = nifti([pathFunc,dirname,filesep,'stmpW_',fileFunc{yy}]);
    finalEPI = preFiReStru.dat(:,:,:,:);
    clear preFiReStru
    delete([pathFunc dirname filesep 'stmpW_' fileFunc{yy}])
    delete([pathFunc dirname filesep 'tmpW_' fileFunc{yy}])
    fprintf('Done! ============== %s\n\r',spm('time'));
    
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
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Creating DVARs Temporal Mask\r\n');
        dvars_TM = zeros(size(finalEPI,4),1);
        dvars_TM = MeanDvars>str2num(get(handles.DVARthres,'String'));
        fprintf('Done! ============== %s\r\n',spm('time'));
        fprintf(Totmovtxt,'%d \t',sum(dvars_TM));
        fprintf(Motfid,'Number of DVAR suprathreshold volumes: %d \r\n',sum(dvars_TM));
        DVAR_FD_BIN(2) = 1;
    else
        fprintf(Totmovtxt,'N/A\t');
        dvars_TM = [];
        DVAR_FD_BIN(2) = 0;
        MeanDvars = [];
        fprintf(Totmovtxt,'N/A\t');
    end

    newEPI = Func;
    newEPI.dat.fname = [pathFunc,dirname,filesep,'ThreSW_',fileFunc{yy}];
    newEPI.descrip = 'UF²C-Masked and Thresholded_sw_EPI';
    newEPI.dat.dtype = 'INT16-LE';
    newEPI.dat(:,:,:,:) = finalEPI;
    create(newEPI)

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
                imgRR = getframe(figu);
                imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Regression Matrix.tif']);
            end
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d regressors creation. Check your data.',yy), 'Process  Aborted')
           fclose('all');
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
    TM_DVARsV = 0;
    TM_FDV = 0;
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
    
    RemVet = FD_TM + dvars_TM;
    RemVet = double(RemVet>0);
    nRVol = sum(RemVet);
    
    fprintf(Motfid,'Total Number of suprathreshold volumes: %d',nRVol);
    fprintf(Totmovtxt,'%d\r\n',nRVol);
    
    fclose(Motfid);
    
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Creating Structure of data\r\n');

    eval(sprintf('StruTOT.Subj%d.Fname_s = [''ThreSW_'',fileFunc{yy}];',yy,yy));
    eval(sprintf('StruTOT.Subj%d.Path_s = [pathFunc,dirname];',yy));
    eval(sprintf('StruTOT.Subj%d.Frgval_s = Frgval;',yy));
    eval(sprintf('StruTOT.Subj%d.Regstr_s = Regstr;',yy));
    eval(sprintf('StruTOT.Subj%d.RemVet_s = RemVet;',yy));
    eval(sprintf('StruTOT.Subj%d.nRVol_s = nRVol;',yy));
    eval(sprintf('StruTOT.Subj%d.FD_TS_s = FDts;',yy));
    eval(sprintf('StruTOT.Subj%d.DVAR_TS_s = MeanDvars;',yy));
    eval(sprintf('StruTOT.Subj%d.DVAR_FD_Bin_s = DVAR_FD_BIN;',yy));

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
    
    fprintf('Done! ============== %s\r\n',spm('time'));

    clear newEPI finalEPIresh Frgval rgval FDts  DVAR_FD_BIN Frgval Regstr RemVet nRVol finalEPI
    
    try
        close(figu)
    end
    try
        close(histCut)
    end

end
fprintf('Imaging Processing and Quantifications: Done! %s\r\n',spm('time'));

multiWaitbar('Total Progress', 0.5);
multiWaitbar('Total Preprocessing Progress',1);

fprintf(Totmovtxt,'\r\n');
fprintf(Totmovtxt,'Volunteers excluded because the number of volumes with suprathresholds movement exceeded the maximum (%s)\r\n',get(handles.LVR,'String'));

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Removing Suprathreshold Volumes\r\n');
fprintf('Step 1:');
StruField = fieldnames(StruTOT);
for yy2 = 1:size(fileFunc,2)
    nOfExclu(yy2) =  eval(sprintf('StruTOT.%s.nRVol_s;',StruField{yy2}));
end
nOfExclu2 = nOfExclu;
nOfExclu2(nOfExclu2>str2num(get(handles.LVR,'String'))) = [];
maxExclu2 = max(nOfExclu2);
fileFunc2 = fileFunc;

fprintf(' Done!\r\n');
fprintf('Step 2: ');

for yy2 = 1:size(fileFunc,2)
    nOfExclus =  eval(sprintf('StruTOT.%s.nRVol_s;',StruField{yy2}));
    Frgval = eval(sprintf('StruTOT.%s.Frgval_s;',StruField{yy2}));
    DVAR_FD_Bin = eval(sprintf('StruTOT.%s.DVAR_FD_Bin_s;',StruField{yy2}));
    
    if DVAR_FD_Bin(2)
        DVARexc = eval(sprintf('StruTOT.%s.DVAR_TS_s;',StruField{yy2}));
        [Btt,Itt] = sort(DVARexc,'descend');
    else DVAR_FD_Bin(1)
        FDexc = eval(sprintf('StruTOT.%s.FD_TS_s;',StruField{yy2}));
        [Btt,Itt] = sort(FDexc,'descend');
    end
    
    mkdir([pathFunc '1-Excluded_Volunteers'])
        
    if nOfExclus<str2num(get(handles.LVR,'String'))
        StruQ = nifti([eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep eval(sprintf('StruTOT.%s.Fname_s',StruField{yy2}))]);
        RemVols = eval(sprintf('StruTOT.%s.RemVet_s',StruField{yy2}));
        RemVols = find(RemVols);
        
        FinalRem = [RemVols' Itt(1:maxExclu2)'];
        
        try
            FinalRem = unique(FinalRem,'legacy');
        catch
            fprintf('Switching to "OLD Matlab version" option\r\n')
            [Cs,ia,tth] = unique(FinalRem,'first');
            ia = sort(ia);
            for iop = 1:numel(Cs)
                FRtmp(iop) = FinalRem(ia(iop));
            end
            FinalRem = FRtmp;
            fprintf('Done!\r\n')
            fprintf('Switching back\r\n')
        end
        
        if numel(FinalRem) > maxExclu2
            FinalRem = FinalRem(1:maxExclu2);
        end
        
        mat101 = StruQ.dat(:,:,:,:);
        mat101(:,:,:,FinalRem) = [];
        if ~isempty(Frgval)
            Frgval(FinalRem,:) = [];
        end
        
        eval(sprintf('StruTOT.%s.Frgval_s = Frgval;',StruField{yy2}));
        
        newEPI = Func;
        funcNAME = eval(sprintf('StruTOT.%s.Fname_s',StruField{yy2}));
        newEPI.dat.fname = [eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep 'VolResThreSW_' funcNAME(8:end)];
        newEPI.descrip = 'UF²C-VolRemoved-Thresholded_sw_EPI';
        newEPI.dat.dtype = 'INT16-LE';
        newEPI.dat.dim = size(mat101);
        newEPI.dat(:,:,:,:) = mat101;
        create(newEPI)
    else
        [pathTrash,dirnameTMP] = fileparts(eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})));
        movefile(eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})),[pathFunc '1-Excluded_Volunteers' filesep dirnameTMP])
        fprintf(Totmovtxt,'%s \r\n',dirnameTMP);
        fileFunc2{yy2} = 'Excluded';
        eval(sprintf('StruTOT.%s.Path_s = ''Excluded'';',StruField{yy2}));
    end
end
fprintf('Done!\r\n');

fprintf('Step 3: ');

idxs = 1;
IncluSubj = [];
for yy2 = 1:size(fileFunc2,2)
    if ~strcmp(fileFunc2{yy2},'Excluded')
        IncluSubj = [IncluSubj yy2];
    end
end

fprintf('Done!\r\n');
fprintf('Done! ============== %s\r\n',spm('time'));

multiWaitbar('Total Progress', 0.75);
multiWaitbar('Total Filtering & Regression Process', 0, 'Color', 'b' ); %CHARGE WAIT BAR

multiWaitbar('Filtering & Regression', 0, 'Color', 'w' ); %CHARGE WAIT BAR

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Performing regressions and/or filtering\r\n');

for yy2 = IncluSubj
    multiWaitbar('Total Filtering & Regression Process', (yy2/size(fileFunc2,2))*0.95 );
    multiWaitbar('Filtering & Regression', 'Reset' );
    
    if size(fileFunc2{yy2},2)<15
        fprintf('Subj. %d  (%s...)\n',yy2,fileFunc2{yy2});
    else
        fprintf('Subj. %d  (%s...)\n',yy2,fileFunc2{yy2}(1:14));
    end

    StruPre = nifti([eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep 'VolResThreSW_' fileFunc2{yy2}]);
    finalEPI = StruPre.dat(:,:,:,:);
    Dt = StruPre.dat.dim(4);
    Regstr = eval(sprintf('StruTOT.%s.Regstr_s',StruField{yy2}));
    Frgval = eval(sprintf('StruTOT.%s.Frgval_s',StruField{yy2}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtering and regressions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Filtering and Regressing === %s\n',spm('time'));
    if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
        BinFazFilt = get(handles.FazFilt,'Value');
        finalEPI = filtregr_uf2c(finalEPI,BinFazFilt,Dx,Dy,Dz,Dt,Frgval,HdLOW,HdHIGH,aFil,tmp1,tmp2,Regstr);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creating the FiltRegre Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Creating and Smoothing the Final Image === %s\n',spm('time'));
    try
        if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
            if  get(handles.FazFilt,'Value') && ~isequal(Regstr,'')
                newEPI = Func;
                newEPI.dat.fname = [eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep,'FiltRegrSW_VRem_',fileFunc2{yy2}];
                newEPI.descrip = 'UF²C-Filtered_Regressed_sw_EPI';
                newEPI.dat.dtype = 'INT16-LE';
                newEPI.dat.dim = size(finalEPI);
                newEPI.dat(:,:,:,:) = finalEPI;
                create(newEPI)
                
                fprintf('Done! ============== %s\n\r',spm('time'));
                
            else
                if get(handles.FazFilt,'Value') && isequal(Regstr,'')
                    newEPI = Func;
                    newEPI.dat.fname = [eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep,'FiltSW_VRem_',fileFunc2{yy2}];
                    newEPI.descrip = 'UF²C-Filtered_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat.dim = size(finalEPI);
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                else
                    newEPI = Func;
                    newEPI.dat.fname = [eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep,'RegrSW_VRem_',fileFunc2{yy2}];
                    newEPI.descrip = 'UF²C-Regressed_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat.dim = size(finalEPI);
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                end
            end
        else
            newEPI = Func;
            newEPI.dat.fname = [eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep,'ThreSW_VRem_',fileFunc2{yy2}];
            newEPI.descrip = 'UF²C-Masked and Thresholded_sw_EPI';
            newEPI.dat.dtype = 'INT16-LE';
            newEPI.dat.dim = size(finalEPI);
            newEPI.dat(:,:,:,:) = finalEPI;
            create(newEPI)
            fprintf('Done! ============== %s\n\r',spm('time'));
        end
    catch
        set(handles.status,'String','Error...')
        warndlg(sprintf('An error occured during subject %d image creation. Check the avaiable disk space.',yy2), 'Process   Aborted')
       fclose('all');
        return
    end
    
    delete([eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep 'ThreSW_' fileFunc2{yy2}])
    delete([eval(sprintf('StruTOT.%s.Path_s',StruField{yy2})) filesep 'VolResThreSW_' fileFunc2{yy2}])
    
    fprintf('Subject %d Done! ====== %s\n\r',yy2,spm('time'));
    
end
multiWaitbar('Total Filtering & Regression Process', 1);
multiWaitbar('Total Progress', 1);
multiWaitbar('CloseAll'); %close the progressbar window

set(handles.status,'String','Done!!!')

if license('test','Distrib_Computing_Toolbox')
    ps.Pool.AutoCreate = PsOri;
end

fclose('all');

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
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)

else
    if get(handles.extWM,'Value') || get(handles.extCSF,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)

    end
    if isequal(get(handles.extCSF,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extGM,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    end
end
    
function extWM_Callback(hObject, eventdata, handles)
if get(handles.extWM,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extCSF,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extCSF,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extWM,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    end

end

function extCSF_Callback(hObject, eventdata, handles)
if get(handles.extCSF,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extWM,'Value')
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)

    end
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.FDcheck,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.extCSF,'Value',1)
        set(handles.text27,'Enable','on')
        set(handles.DVARthres,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)

    end
end

function FDcheck_Callback(hObject, eventdata, handles)
if get(handles.FDcheck,'Value')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
    if isequal(get(handles.extGM,'Value'),0) && isequal(get(handles.extCSF,'Value'),0) && isequal(get(handles.extWM,'Value'),0)
        warndlg('You should to select at leat one motion quantification parameter','Attention')
        set(handles.FDcheck,'Value',1)
        set(handles.text28,'Enable','on')
        set(handles.FDthres,'Enable','on')
    end
end

function LVR_Callback(hObject, eventdata, handles)

function LVR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DVAR_TS_Callback(hObject, eventdata, handles)

function EPCo_Callback(hObject, eventdata, handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
global PsOri ps
if license('test','Distrib_Computing_Toolbox')
    ps.Pool.AutoCreate = PsOri;
end
delete(hObject);
