function varargout = FuncConSW(varargin)
% UF²C M-file for FuncInte.fig
% UF²C - User Friendly Functional Connectivity
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
                   'gui_OpeningFcn', @FuncConSW_OpeningFcn, ...
                   'gui_OutputFcn',  @FuncConSW_OutputFcn, ...
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

function FuncConSW_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = FuncConSW_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

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
fprintf('==========================================\r\n');

roicoord = [str2num(get(handles.ROIinp,'String'))]; % Center ROI coordinate (MNI152) (PCC)

seedSizeMat = [str2num(get(handles.SSInp,'String'))];

x = seedSizeMat(1); % ROI size in voxels (HALF VALUE FOR EACH SIDE)
y = seedSizeMat(2);
z = seedSizeMat(3);

sWS = str2num(get(handles.wsIn,'String')); %Window Size (SIZE OF THE MOVING AVERAGE TO CALC DE CORRELATIONS)

TR = str2num(get(handles.TRInp,'String'));    % repetition time (s)

if isequal(get(handles.FazPrep,'Value'),0)
    if ~isequal(size(fileStru,2),size(fileFunc,2))
        warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
          'You need to add one functional image for each strutural and vice versa.';...
          'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
          'Attention: process aborted!');
        return
    end
end

if sWS<16
    respQ = questdlg('The "Sliding Window Size" is smaller than the canonical HRF duration. Do you wat to continue?','Attention!','Yes','Abort','Abort');
    switch respQ
        case 'Yes'
            
        case 'Abort'
            set(handles.status,'String','Process Aborted')
            return
    end
end

multiWaitbar('Total Progress', 0, 'Color', 'g' );  %CHARGE WAIT BAR
if isequal(get(handles.FazPrep,'Value'),0)
    multiWaitbar('Filtering & Regression', 0, 'Color', 'y' ); %CHARGE WAIT BAR
end
multiWaitbar('Sliding Window', 0, 'Color', 'b' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('1-Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

fideSJ = fopen([pathFunc,filesep,foldername,filesep,'1-Subjects.txt'],'w+'); 
for yTy = 1:size(fileFunc,2)
    fprintf(fideSJ,'%d: \t%s\r\n',yTy,fileFunc{1,yTy});
end
fclose(fideSJ);

imgRR = getframe(FuncConSW);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

fideG = fopen([pathFunc,filesep,foldername,filesep,'Total_Results.txt'],'w+'); % CHARGE OUTPUT LOG FILE
    
fprintf(fideG,'                                          Functional Conectivity ROI Analysis\r\n\r\n');
fprintf(fideG,'Number of subjects included: %d \r\n\r\n', size(fileFunc,2));
fprintf(fideG, 'Subject_file_name');
fprintf(fideG, '                                        Average_Positive_Correlation_Value       Average_Negative_Correlation_Value      Average_Positive_Right_Value      Average_Positive_Left_Value      Average_Negative_Right_Value      Average_Negative_Left_Value\r\n');

if isequal(get(handles.FazPrep,'Value'),0) % just if you need preprocessing
    if isequal(get(handles.FazFilt,'Value'),1)

        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Designing filters\n');
        %%%%%%%%%%%%%%%%%%%%%% FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        fprintf('Done! ============== %s\r\n',spm('time'));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WinAvSavePos = struct;
WinAvSaveNeg = struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF SUBJECT LOOP %%%%%%%%%%%%%%%%%%%

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);
SPMdir2 = which('spm');
SPMdir2 = SPMdir2(1:end-5);

for yy = 1: size(fileFunc,2)
    
    multiWaitbar('Total Progress', (yy/size(fileFunc,2))*0.95 );
    
    fprintf('\r\n')
    fprintf('UF²C ================= %s\n',spm('time'));
    fprintf('Starting subj. %d  (%s...) processing\n',yy,fileFunc{yy}(1:10));
    fprintf('============================================\r\n\r\n');
    
    try
        close(fig)
    end
    
    Funpre = nifti([pathFunc,fileFunc{yy}]);
    
    fvs = double(Funpre.hdr.pixdim(2:4)); %GET PIXEL SIZE (FUNCTIONAL)
    
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
        for tt = 1:Funpre.dat.dim(4)
            file1{tt} = [pathFunc,dirname,filesep,fileFunc{yy},',',sprintf('%d',tt)];
        end
    else
        multiWaitbar('Sliding Window', 'Reset' );
        dirname = Funpre.dat.fname;
        [ax,dirname,cx] = fileparts(dirname);
        dirname = dirname(12:end);
        
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

        mkdir(pathFunc,dirname)
    end
    
    mkdir([pathFunc,dirname],'Correlation_map')
    mkdir([pathFunc,dirname,filesep,'ROI_mask'])
    mkdir([pathFunc,dirname],'Results_Log')
    
    fide = fopen([pathFunc,dirname,filesep,'Results_Log',filesep,'Results_Log.txt'],'w+');
    
    if isequal(get(handles.FazPrep,'Value'),0)
        try
             BounB = get(handles.SPMbb,'Value');
             Preproc_uf2c(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname)
        catch
             set(handles.status,'String','Error...')
             warndlg(sprintf('An error occured during subject %d preprocessing. This is a SPM error. Check your data',yy), 'Proces Aborted')
             return
        end


        try
            close(fig);
        end
        
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Extracting fMRI global signals\n');
        
        Func = nifti([pathFunc,dirname,filesep,'sw',fileFunc{yy}]);
        Struct = nifti([pathFunc,dirname,filesep,'wm',fileStru{yy}]);

        matrix4 = Func.dat(:,:,:,:);
        matrix3 = Struct.dat(:,:,:);

        sizeFunc = size(matrix4);

        template = matrix4;
        template(:,:,:) = 1;
        templ = Func;
        templ.dat.dim = [sizeFunc(1) sizeFunc(2) sizeFunc(3)];
        templ.dat.fname = [pathFunc,dirname,filesep,'Template.nii'];
        templ.dat.dtype = 'FLOAT32-LE';
        templ.dat(:,:,:) = template;
        create(templ)

        %%%%%%%%%%%%%%%%%%%%%%% Interp segmented images %%%%%%%%%%%%%%%%%%%%%%%%%%
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %GM image
                                            [pathFunc,dirname,filesep,'wc1',fileStru{yy}]};
        matlabbatch{1}.spm.util.imcalc.output = 'c1interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.15';% Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)
        
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc,dirname,filesep,'c1interp.nii']};
        if get(handles.SPMbb,'Value')
            matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        else
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
        end
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch)
        
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'sc1interp.nii']};
        matlabbatch{1}.spm.util.imcalc.output = 'c1interpF.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.3';% Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)

        c1Nii = nifti([pathFunc,dirname,filesep,'c1interpF.nii']);
        c1mat = c1Nii.dat(:,:,:);

        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %WM image
                                            [pathFunc,dirname,filesep,'wc2',fileStru{yy}]};
        matlabbatch{1}.spm.util.imcalc.output = 'c2interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.3';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)

        c2Nii = nifti([pathFunc,dirname,filesep,'c2interp.nii']);
        c2mat = c2Nii.dat(:,:,:);

        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %CSF image
                                                [pathFunc,dirname,filesep,'wc3',fileStru{yy}]};
        matlabbatch{1}.spm.util.imcalc.output = 'c3interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.3';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)
        c3Nii = nifti([pathFunc,dirname,filesep,'c3interp.nii']);
        c3mat = c3Nii.dat(:,:,:);

        MeanWM = zeros(1,sizeFunc(4));   % Prelocation
        MeanCSF = zeros(1,sizeFunc(4)); % Prelocation
        finalEPI = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c2Map = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c3Map = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation

        for vx = 1:sizeFunc(4); 
            finalEPI(:,:,:,vx) = c1mat.*matrix4(:,:,:,vx);% Apply Grey Matter mask to the original EPI to remove outside brain voxels

            c2Map(:,:,:,vx) = c2mat.*matrix4(:,:,:,vx); % Apply White Matter mask to the original EPI. Regress to WM variability
            dataWM = c2Map(:,:,:,vx);
            NumeWM = sum(sum(sum(dataWM)));  
            DenWM = nnz(dataWM);
            MeanWM(vx) = NumeWM/DenWM; %Calculates the average WM time series.

            c3Map(:,:,:,vx) = c3mat.*matrix4(:,:,:,vx); % Apply CSF mask to the original EPI. regress to CSF variability
            dataCSF = c3Map(:,:,:,vx);
            NumeCSF = sum(sum(sum(dataCSF)));  
            DenCSF = nnz(dataCSF);
            MeanCSF(vx) = NumeCSF/DenCSF; %Calculates the average CSF time series.
        end
        
        clear c1mat c2mat c3mat c1Nii c2Nii c3Nii
        fprintf('Done! ============== %s\r\n',spm('time'));

        %nvol23 = str2double(get(handles.numDyna,'String'));
        r1 = '';
        r2 = '';
        r3 = '';

        if isequal(get(handles.checkMOV,'Value'),1)
            rpfile = textread([pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt']); % read alignment file output
            rgval1 = rpfile; %add movement parameters to regress matrix

            Min_rp_1 = min(rgval1(2:end,1));
            Min_rp_2 = min(rgval1(2:end,2));
            Min_rp_3 = min(rgval1(2:end,3));
            Min_rp_4 = min(rgval1(2:end,4));
            Min_rp_5 = min(rgval1(2:end,5));
            Min_rp_6 = min(rgval1(2:end,6));

            rgval1(1,:) = 0;

            rgval1(2:end,1) = rgval1(2:end,1)+ abs(Min_rp_1);
            rgval1(2:end,2) = rgval1(2:end,2)+ abs(Min_rp_2);
            rgval1(2:end,3) = rgval1(2:end,3)+ abs(Min_rp_3);
            rgval1(2:end,4) = rgval1(2:end,4)+ abs(Min_rp_4);
            rgval1(2:end,5) = rgval1(2:end,5)+ abs(Min_rp_5);
            rgval1(2:end,6) = rgval1(2:end,6)+ abs(Min_rp_6);

            r1 = 'mov-';
        end

        if isequal(get(handles.checkOsc,'Value'),1)
            rgval2 = transpose([MeanWM;MeanCSF]); %add mean white matter and CSF parameter to regress matrix
            r2 = 'osc-';
        end

        if isequal(get(handles.checkReg,'Value'),1)
            if ~isequal(nrg,0)
                rgval3 = load([tmpDIR 'additional_Reg.mat']);
                nCol = 1;
                for yt = 1:numofVoluF:nrg*numofVoluF

                    rgval3r(:,nCol) = rgval3.regVets(yt:yt+(numofVoluF-1),yy);

                    nCol = nCol +1;
                end
                r3 = 'addt';

            end
        end

        if isequal(r3,'addt')
            if size(unique(rgval3r),1)==1
                r3 = '';
            end
        end

        Regstr = [r1,r2,r3];

        if isequal(Regstr,'osc-')
            rgval = rgval2;
        end
        if isequal(Regstr,'mov-')
            rgval = rgval2;
        end
        if isequal(Regstr,'addt')
            rgval = rgval2;
        end
        if isequal(Regstr,'mov-osc-')
            rgval = [rgval1,rgval2];
        end
        if isequal(Regstr,'mov-addt')
            rgval = [rgval1,rgval3r];
        end
        if isequal(Regstr,'mov-osc-addt')
            rgval = [rgval1,rgval2,rgval3r];
        end
        if isequal(Regstr,'osc-addt')
            rgval = [rgval2,rgval3r];
        end
        if isequal(Regstr,'')
            warndlg('The use of regressos are strongly recommended!','Attention');
        end

        SignCorr = zeros(size(finalEPI,4),1); %Prelocation

        if ~isequal(Regstr,'')
            for rte = 1:size(rgval,2) % all regressors are normalized: minimal value = 0 maximal value = 1. the fiorst line is a copy of the second.
                MinV = min(rgval(2:end,rte));
                MaxV = max(rgval(2:end,rte));
                for ttt = 2:size(finalEPI,4)
                    Frgval(ttt,rte) = ((rgval(ttt,rte))-MinV)/(MaxV-MinV);
                end
            end

            Frgval = [Frgval,ones(sizeFunc(4),1)];
            Frgval(1,:) = Frgval(2,:);

            figu = figure('Visible','on');
            set(figu,'Name','RegressionMatrix',...
                'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.3 ScreSize(1)*.4]),...
                'Color',[1 0.94 0.86]);
            imagesc(Frgval);
            grid('off')

            colormap(gray)
            title('Regression Matrix','FontSize', 14);
            drawnow
        %             saveas(figu,[pathFunc,dirname,filesep,'Correlation_map',filesep,'Regression Matrix'],'tif')
            imgRR = getframe(gcf);
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Correlation_map',filesep,'Regression Matrix.tif']);
            fprintf('Done! ============== %s\r\n',spm('time'));    
        end

        Dx = size(finalEPI,1);
        Dy = size(finalEPI,2);
        Dz = size(finalEPI,3);
        Dt = size(finalEPI,4);

        THREimg = mean(finalEPI,4);
        THREimgresh = reshape(THREimg,prod([Dx,Dy]),Dz)';
        THREimgresh(THREimgresh==0) = NaN;

        try %condition to old Matlab versions
            minVV = min(min(THREimgresh));
            maxVV = max(max(THREimgresh));
            nbins = round((maxVV-minVV)/5);
            [Nxs,edges,binxs] = histcounts(THREimgresh,nbins);
            [Mxs,Ixs] = max(Nxs);
            HigherFreq = edges(Ixs);
            CutoffValue = round(0.33*HigherFreq); %% Threshold: The cutoff intensity is 33% of the most frequent intensity
        catch
            disp('--- Old Matlab Version: Alternative code\r\n')
            THREimgresh = reshape(THREimgresh,prod([Dx,Dy,Dz]),1)';
            minVV = min(THREimgresh);
            maxVV = max(THREimgresh);
            nbins = ceil((maxVV-minVV)/5);
            [Nxs,edges] = hist(THREimgresh,nbins);
            [Mxs,Ixs] = max(Nxs);
            HigherFreq = edges(Ixs);
            CutoffValue = round(HigherFreq*0.33);
        end

        histCut = figure('Visible','on');
        set(histCut,'Name','Average Image Histogram',...
            'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.4 ScreSize(1)*.3]),...
                'Color',[1 0.94 0.86]);
        try %condition to old Matlab versions
            histogram(THREimgresh,nbins);
        catch
            bar((minVV:5:maxVV),Nxs,'hist');
        end

        xlabel('Intensity'  ,'FontSize',14);
        ylabel('Frequency','FontSize',14);
        hold on
        ylimI = get(gca,'YLim');
        line([CutoffValue,CutoffValue],ylimI,'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
        str1 = ('\leftarrow');
        str2 = sprintf('Cutoff Intensity:\n %d',CutoffValue);
        strFi = [str1 str2];
        text(CutoffValue,ylimI(2)/2,strFi)
        title('Average Image Histogram','FontSize', 14);
        drawnow
    %         saveas(histCutoff,[pathFunc,dirname,filesep,'Correlation_map',filesep,'Avg_Img_Histogram'],'tif')
        imgRR = getframe(gcf);
        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Correlation_map',filesep,'Avg_Img_Histogram.tif']);

        clear minVV maxVV str1 str2 strFi histCutoff THREimgresh Nxs edges THREimg

        finalEPI(finalEPI<CutoffValue) = 0;

        if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')

            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Performing regressions and/or filtering\r\n');

            finalEPIresh = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';
            vetBin = (sum(finalEPIresh)~=0);
            finalEPIresh(:,~any(finalEPIresh,1) ) = [];  %columns

            NumT1 = Dt;% number of time points
            aFil = 1;

            for tt3 = 1:size(finalEPIresh,2)
                multiWaitbar('Filtering & Regression', tt3/size(finalEPIresh,2));
                SignCorr = zeros(Dt,1);
                SignCorr = finalEPIresh(:,tt3);
                if ~isequal(Regstr,'')
                    b = regress(SignCorr,Frgval); %Regression. Includes: Movement regression(6 parameters), mean White Matter signal and mean CSF signal
                    s_new  = Frgval*b;
                    s_corr = SignCorr - s_new;
                else
                    s_corr = SignCorr;
                end
                if isequal(get(handles.FazFilt,'Value'),1)
                    
                    tmp1 = 3*(max(length(HdLOW.Numerator),length(aFil))-1);
                    tmp2 = 3*(max(length(HdHIGH.Numerator),length(aFil))-1);

                    %Low Pass Filtering
                    if NumT1<=tmp1
                        s_filt2 = repmat(s_corr,ceil(tmp1/NumT1),1);
                    else
                        s_filt2 = s_corr;
                    end
                    s_filtAB1LOW = filtfilt(HdLOW.Numerator,aFil,s_filt2);
                    s_filtAB1LOW = s_filtAB1LOW(1:NumT1);

                    %High Pass Filtering
                    if NumT1<=tmp2
                        s_filtAB1LOW = repmat(s_filtAB1LOW,ceil(tmp2/NumT1),1);
                    end
                    s_filtAB1LOW_HIGH = filtfilt(HdHIGH.Numerator,aFil,s_filtAB1LOW);
                    s_filtAB1LOW_HIGH = s_filtAB1LOW_HIGH(1:NumT1);

                    %Detrend
                    s_filtABF = detrend(s_filtAB1LOW_HIGH,'linear');  %detrend with breakpoints over 128 seconds
                    finalEPIresh(:,tt3) = s_filtABF(:);
                else
                    finalEPIresh(:,tt3) = s_corr;
                end
            end
            clear s_filtAB1LOW s_filtAB1LOW_HIGH s_filtABF s_filt2 s_corr s_new file1 matrix4 c2Map c3Map
            finalEPIresh2 = zeros(Dt,size(vetBin,2));

            finalEPIresh2 = bsxfun(@and,vetBin,ones(1,Dt)');
            finalEPIresh2 = double(finalEPIresh2);
            finalEPIresh2(finalEPIresh2~=0)=finalEPIresh;

            finalEPI = reshape(finalEPIresh2',Dx,Dy,Dz,Dt);
            finalEPI = squeeze(finalEPI);
            finalEPI = round(1000.*finalEPI);
            fprintf('Done! ============== %s\r\n',spm('time'));
        end
        newEPI = Func;
        newEPI.dat.fname = [pathFunc,dirname,filesep,'FiltRegrSW_',fileFunc{yy}];
        newEPI.descrip = 'UF²C-Filtered_Regressed_sw_EPI';
        newEPI.dat.dtype = 'INT16-LE';
        newEPI.dat(:,:,:,:) = finalEPI;
        create(newEPI)
        clear newEPI finalEPIresh Frgval rgval finalEPIresh2

        if exist('figu')
            close(figu)
        end
        if exist('histCut')
            close(histCut)
        end
    else
        Func = nifti([pathFunc,fileFunc{yy}]);
        finalEPI = Func.dat(:,:,:,:);
        Dx = size(finalEPI,1);
        Dy = size(finalEPI,2);
        Dz = size(finalEPI,3);
        Dt = size(finalEPI,4);
    end

    clear('finalEPIresh')
    clear('Frgval')
    clear('rgval')
    
%   Convert from MNI coord to matrix coord.
    if get(handles.FazPrep,'Value')
        Vhdr = spm_vol([pathFunc,filesep,fileFunc{yy}]);
    else
        Vhdr = spm_vol([pathFunc,dirname,filesep,'FiltRegrSW_',fileFunc{yy}]);
    end
    
    TcM = Vhdr.mat;
    EPIcoord = [roicoord(:,1) roicoord(:,2) roicoord(:,3) ones(size(roicoord,1),1)]*(inv(TcM))';
    EPIcoord(:,4) = [];
    EPIcoord = round(EPIcoord);
    
    Sxx = x/2;
    Syy = y/2;
    Szz = z/2;
    
    if Sxx<1
        Sxx=0;
    end
    if Syy<1
        Syy=0;
    end
    if Szz<1
        Szz=0;
    end
    
    roi1 = zeros(Dx,Dy,Dz);%%%%%%%%%%%%%%%%%%%% Create A cubic ROI centered in the input coordinate %%%%%%%%%%%%%%%%%%%%%%
    roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3):EPIcoord(3)+Szz) = 1;
    roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3):EPIcoord(3)+Szz) = 1;
    
    roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3):EPIcoord(3)+Szz) = 1;
    roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3)-Szz:EPIcoord(3)) = 1;
    
    roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3)-Szz:EPIcoord(3)) = 1;
    roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3)-Szz:EPIcoord(3)) = 1;
    
    roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3)-Szz:EPIcoord(3)) = 1;
    roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3):EPIcoord(3)+Szz) = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [iCo,jCo,kCo]= ind2sub(size(roi1), find(roi1>0));
    coords = [iCo,jCo,kCo];
    
    for vx = 1:Dt; % get the ROI voxels time series. This step exclude from the ROI eventuals white matter voxels.
        final3D(:,:,:) = roi1.*finalEPI(:,:,:,vx);
        for hj = 1:size(coords,1)
            vetsCor(vx,hj) = final3D(coords(hj,1),coords(hj,2),coords(hj,3));
        end
    end
    
    verifCor = mean(vetsCor,2);
    xxX = corr(vetsCor(:,:),verifCor);
    xxX(isnan(xxX)) = 0;

    xxX_Z = 0.5.*(log(1+xxX) - log(1-xxX)); %Z transf
    stdxxX = std(xxX_Z);
    meanxxX = mean(xxX_Z);
    lowerCut = meanxxX-stdxxX;
    alower = (xxX_Z<lowerCut);  

    CutFinal = alower;
    CutFinal = (abs(CutFinal-1))';
    CutFinal = repmat(CutFinal,size(vetsCor,1),1);
    vetsCorF = CutFinal.*vetsCor;
    vetsCorF(:,~any(vetsCorF,1) ) = [];
    Size4D = size(vetsCorF,2);
    mean4D = (mean(vetsCorF,2))';
    fprintf(fide,'%d voxel(s) was (were) removed from the \r\n', (size(coords,1)-Size4D));
    fprintf(fide,'seed due hight temporal variability when compared to the anothers \r\n\r\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Correlation %%%%%%%%%%%%%%%%%%%%%
    finalEPI = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';
    map = zeros(Dt-sWS,prod([Dx,Dy,Dz]));%Prelocation
    Vn = 1;
    %sWS = 16;
    
    
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Estimating correlations\n');

    for MA = 1:Dt-sWS
        multiWaitbar('Sliding Window', MA/(Dt-sWS) );
        map(Vn,:) = corr(finalEPI(MA:MA+(sWS-1),:),mean4D(MA:MA+(sWS-1))');
        Vn = Vn+1;
    end
    map(isnan(map)) = 0;
    map = reshape(map',Dx,Dy,Dz,Dt-sWS);
    map = squeeze(map);
    fprintf('Done! ============== %s\r\n',spm('time'));
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Estimating averages\n');

    mapP = map;
    mapP(mapP(:,:,:,:)<0)=0;
    mapN = map;
    mapN(mapN(:,:,:,:)>0)=0;
    
    Mean4DVetPOS = zeros(1,(Vn-1));
    Mean4DVetNEG = zeros(1,(Vn-1));
    
    for ix = 1:(Vn-1)
        baseP = mapP(:,:,:,ix);
        NumCorWindowP = sum(sum(sum(baseP)));
        DenCorWindowP = nnz(baseP);
        corrWindowP = NumCorWindowP/DenCorWindowP;
        Mean4DVetPOS(ix) = corrWindowP;
        
        baseN = mapN(:,:,:,ix);
        NumCorWindowN = sum(sum(sum(baseN)));
        DenCorWindowN = nnz(baseN);
        corrWindowN = NumCorWindowN/DenCorWindowN;
        Mean4DVetNEG(ix) = corrWindowN;
    end
    
    meanPosR = mean(Mean4DVetPOS);
    meanNegR = mean(Mean4DVetNEG);
    fprintf('Done! ============== %s\r\n',spm('time'));

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Plotting\n');

    fig = figure;
    set(fig,'Name',dirname,...
        'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1).*0.6 ScreSize(1).*0.4]),...
                    'Color',[1 0.94 0.86]);
    subplot1 = subplot(1,2,1);
    plot(Mean4DVetPOS,'Parent',...
                    subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point' ,'FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Positive Moving Window Average Cor value');
    
    subplot2 = subplot(1,2,2);
    plot(Mean4DVetNEG,'Parent',...
                    subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    
    title('Negative Moving Window Average Cor value');
    
    Y1 = Mean4DVetPOS;
    Y2 = Mean4DVetNEG;

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
    title('Total Negative Moving Window Average Corrrelation value');
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
    
    fprintf('Done! ============== %s\r\n',spm('time'));

    eval(sprintf('WinAvSavePos.suj_%d = Mean4DVetPOS',yy));
    eval(sprintf('WinAvSaveNeg.suj_%d = Mean4DVetNEG',yy));
        
    save([pathFunc,dirname,filesep,'Correlation_map',filesep,'Positive_Serie'],'Mean4DVetPOS')
    save([pathFunc,dirname,filesep,'Correlation_map',filesep,'Negative_Serie'],'Mean4DVetNEG')

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Creating images\n');

    RoiNii = Func;   % Creates a NIfTI object
    WinCorr4DP = Func; 
    WinCorr4DN = Func; 
    
    RoiNii.dat.dim = [Dx Dy Dz]; % Apply the correct matrix size
    WinCorr4DP.dat.dim = [Dx Dy Dz (Vn-1)];
    WinCorr4DN.dat.dim = [Dx Dy Dz (Vn-1)];

    RoiNii.dat.fname = [pathFunc,dirname,filesep,'ROI_mask',filesep,dirname,'_Roi.nii']; % Change file name
    WinCorr4DP.dat.fname = [pathFunc,dirname,filesep,'Correlation_map',filesep,dirname,'_Cr4D_P.nii'];
    WinCorr4DN.dat.fname = [pathFunc,dirname,filesep,'Correlation_map',filesep,dirname,'_Cr4D_N.nii'];
    
    RoiNii.dat.dtype = 'FLOAT32-LE'; % verify matrix datatype
    WinCorr4DP.dat.dtype = 'FLOAT32-LE';
    WinCorr4DN.dat.dtype = 'FLOAT32-LE';
    
    RoiNii.descrip = dirname; % Apply Name information in the Nii struct (description field)
    WinCorr4DP.descrip = dirname;
    WinCorr4DN.descrip = dirname;

    RoiNii.dat(:,:,:) = roi1; % add matriz value to the NIfTI object
    WinCorr4DP.dat(:,:,:,:) = mapP;
    WinCorr4DN.dat(:,:,:,:) = mapN;
    
    create(RoiNii) % creates the NIfTI file from the Object
    create(WinCorr4DP)
    create(WinCorr4DN)
    
    fprintf('Done! ============== %s\r\n',spm('time'));

    clear matlabbatch %Start to creates the smoothed correlation 4D map
    
    Corr4D_Mult = cell((Vn-1),1);
    
    for tx=1:(Vn-1)
        Corr4D_MultP{tx,1} = [pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_P.nii' ',' num2str(tx)];
    end
    
    matlabbatch{1}.spm.spatial.smooth.data = Corr4D_MultP;                                       
    if get(handles.SPMbb,'Value')
        matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
    else
        matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    end
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
    spm_jobman('run',matlabbatch)

    
    clear matlabbatch %Start to creates the smoothed correlation 4D map
    Corr4D_Mult = cell((Vn-1),1);
    
    for tx=1:(Vn-1)
        Corr4D_MultN{tx,1} = [pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_N.nii' ',' num2str(tx)];
    end
    
    matlabbatch{1}.spm.spatial.smooth.data = Corr4D_MultN;                                       
    if get(handles.SPMbb,'Value')
        matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
    else
        matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    end
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
    spm_jobman('run',matlabbatch)
    
    equa = 'i1';  %creates the equation to use in ImCalc
    if Dt>Dt-sWS
        for eqn = 2:Dt-sWS
            equa = sprintf('%s+i%d',equa,eqn);
        end
        equa = ['(' equa ')' sprintf('./%d',eqn)];
    
        Corr4D_MultP = transpose(Corr4D_MultP);
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = Corr4D_MultP;
        matlabbatch{1}.spm.util.imcalc.output = [dirname,'_MeanP.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc dirname filesep 'Correlation_map' filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = equa;
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        spm_jobman('run',matlabbatch)

        Corr4D_MultN = transpose(Corr4D_MultN);
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = Corr4D_MultN;
        matlabbatch{1}.spm.util.imcalc.output = [dirname,'_MeanN.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc dirname filesep 'Correlation_map' filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = equa;
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 1;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc dirname filesep 'Correlation_map' filesep dirname '_MeanP.nii']};                                       
        if get(handles.SPMbb,'Value')
            matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        else
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
        end
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc dirname filesep 'Correlation_map' filesep dirname '_MeanN.nii']};                                       
        if get(handles.SPMbb,'Value')
            matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        else
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
        end
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
        spm_jobman('run',matlabbatch)

        % positive values
            StruPF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_P.nii']);
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
            StruNF = nifti([pathFunc dirname filesep 'Correlation_map' filesep dirname '_Cr4D_N.nii']);   
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
            NumRLP = sum(sum(sum(matlatL(matlatL(:,:,:)>0))));
            DenRLP = nnz(matlatL(:,:,:)>0);
            leftR(rt) = NumRLP/DenRLP;
            %Left values
            matlatR = matlat(:,:,:,rt).*MaskR;
            NumRRP = sum(sum(sum(matlatR(matlatR(:,:,:)>0))));
            DenRRP = nnz(matlatR(:,:,:)>0);
            rightR(rt) = NumRRP/DenRRP;

        % negative values
            %Right values
            matlatL2 = matlat2(:,:,:,rt).*MaskL2;
            NumRLP2 = sum(sum(sum(matlatL2(matlatL2(:,:,:)<0))));
            DenRLP2 = nnz(matlatL2(:,:,:)<0);
            leftR2(rt) = NumRLP2/DenRLP2;
            %Left values
            matlatR2 = matlat2(:,:,:,rt).*MaskR2;
            NumRRP2 = sum(sum(sum(matlatR2(matlatR2(:,:,:)<0))));
            DenRRP2 = nnz(matlatR2(:,:,:)<0);
            rightR2(rt) = NumRRP2/DenRRP2;
        end
    
    end    %in case that there is no 4D volume (no window corr)
        
    if isequal(get(handles.OutputFil, 'Value'),1)
        delete([pathFunc dirname filesep fileFunc{yy}])
        delete([pathFunc dirname filesep fileStru{yy}])
        delete([pathFunc dirname filesep fileFunc{yy}(1:end-4) '.mat'])
        delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_inv_sn.mat'])
        delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_sn.mat'])
        delete([pathFunc dirname filesep 'mean' fileFunc{yy}])
        delete([pathFunc dirname filesep 'c1' fileStru{yy}])
        delete([pathFunc dirname filesep 'c2' fileStru{yy}])
        delete([pathFunc dirname filesep 'c3' fileStru{yy}])
        delete([pathFunc dirname filesep 'm' fileStru{yy}])
        delete([pathFunc dirname filesep 'mwc1' fileStru{yy}])
        delete([pathFunc dirname filesep 'mwc2' fileStru{yy}])
        delete([pathFunc dirname filesep 'mwc3' fileStru{yy}])
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

    fprintf(fide,'%s \r\n\r\n',dirname);   % Write data in the output log file
    fprintf(fide,'Average Positive Correlation Value: %.5f\r\n\r\n',meanPosR);% Write data in the output log file
    fprintf(fide,'Average Negative Correlation Value: %.5f\r\n\r\n',meanNegR);% Write data in the output log file
    fprintf(fide,'Right Side Positive Correlation Value: %.5f\r\n\r\n',rightAverageP);% Write data in the output log file
    fprintf(fide,'Left Side Positive Correlation Value: %.5f\r\n\r\n',leftAverageP);% Write data in the output log file
    fprintf(fide,'Right Side Negative Correlation Value: %.5f\r\n\r\n',rightAverageN);% Write data in the output log file
    fprintf(fide,'Left Side Negative Correlation Value: %.5f\r\n\r\n',leftAverageN);% Write data in the output log file

    fprintf(fideG,'%s:',dirname);% Write data in the output log file
    fprintf(fideG,'                           %.5f                            %.5f                           %.5f                            %.5f                           %.5f                            %.5f\r\n',meanPosR,meanNegR,rightAverageP,leftAverageP,rightAverageN,leftAverageN);% Write data in the output log file
    fclose(fide); %close the writable txt file
    
    clear map rightAverageP leftAverageN leftR rightR meanPosR meanNegR Corr4D_MultN Corr4D_MultP
end 

if size(fileFunc,2)>1
    for iy1 = 1:(Vn-1)
        for ix3 = 1:size(fileFunc,2)
            eval(sprintf('STDPos_%d(%d) = WinAvSavePos.suj_%d(%d);',iy1,ix3,ix3,iy1));
            eval(sprintf('STDNeg_%d(%d) = WinAvSaveNeg.suj_%d(%d);',iy1,ix3,ix3,iy1));
        end
        StdFinalPos(iy1) = std(eval(sprintf('STDPos_%d',iy1)));
        StdFinalNeg(iy1) = std(eval(sprintf('STDNeg_%d',iy1)));
    end
    
    StdFinalPos = StdFinalPos./sqrt(size(fileFunc,2));
    StdFinalNeg = StdFinalNeg./sqrt(size(fileFunc,2));

end


graphfinalPos = WinAvSavePos.suj_1;
graphfinalNeg = WinAvSaveNeg.suj_1;
All_series_Pos(1,:) = WinAvSavePos.suj_1;
All_series_Neg(1,:) = WinAvSaveNeg.suj_1;

if size(fileFunc,2)>1
    for ix2 = 2:size(fileFunc,2)
        eval(sprintf('graphfinalPos = graphfinalPos + WinAvSavePos.suj_%d;',ix2));
        eval(sprintf('graphfinalNeg = graphfinalNeg + WinAvSaveNeg.suj_%d;',ix2));
        eval(sprintf('All_series_Pos(%d,:) = WinAvSavePos.suj_%d;',ix2,ix2));
        eval(sprintf('All_series_Neg(%d,:) = WinAvSaveNeg.suj_%d;',ix2,ix2));
    end
    
    graphfinalPos = graphfinalPos./size(fileFunc,2);
    graphfinalNeg = graphfinalNeg./size(fileFunc,2);
    
    save([pathFunc,filesep,foldername,filesep,'Pos_Average_Correlation_Serie'],'graphfinalPos')
    save([pathFunc,filesep,foldername,filesep,'Neg_Average_Correlation_Serie'],'graphfinalNeg')
    
    save([pathFunc,filesep,foldername,filesep,'All_Pos_Correlation_Series'],'All_series_Pos')
    save([pathFunc,filesep,foldername,filesep,'All_Neg_Correlation_Series'],'All_series_Neg')

end

try
    fig2 = figure;
    set(fig2,'Name','Total Moving Window Average',...
        'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1).*0.6 ScreSize(1).*0.4]),...
                    'Color',[1 0.94 0.86]);
    
    subplot1 = subplot(1,2,1);
    plot(graphfinalPos,'Parent',...
                    subplot1,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point' ,'FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Positive Moving Window Average Corrrelation value');

    subplot2 = subplot(1,2,2);
    plot(graphfinalNeg,'Parent',...
                    subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',1);
    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Negative Moving Window Average Corrrelation value');            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1 = graphfinalPos;
    Y2 = graphfinalNeg;

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

    plot2 = errorbar(Y2,StdFinalNeg/2,'Parent',subplot2,'MarkerEdgeColor',[1 0 0],'MarkerSize',8,...
        'Marker','.','LineWidth',2,'Color',[0 0 1],'DisplayName','data 1');

    xlabel('Window Start point','FontSize',12);
    ylabel('Correlation value','FontSize',12);
    title('Total Negative Moving Window Average Corrrelation value');
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
    drawnow
    imgRR = getframe(gcf);
    imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'Totals_Windows_Graph.tif']);
catch
    
end
    
fclose(fideG); %close the writable txt file
multiWaitbar('CloseAll'); %close the progressbar window
fclose('all');
set(handles.status,'String','Done!!!')
fprintf('\r\n')
fprintf('All Done! ========== %s\n',spm('time'));
fprintf('==========================================\r\n');


function addreg_Callback(hObject, eventdata, handles)
global regVets

try
    clear(regVets)
end
Addreg

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

function SPMbb_Callback(hObject, eventdata, handles)
if isequal(get(handles.SPMbb,'Value'),1)
    set(handles.SSInp,'String','2 2 2')
else
    set(handles.SSInp,'String','4 4 4')
end

function FazPrep_Callback(hObject, eventdata, handles)
if get(handles.FazPrep,'Value')
    set(handles.addStru,'Enable','off')
    set(handles.FazFilt,'Enable','off')
    set(handles.LPF,'Enable','off')
    set(handles.HPF,'Enable','off')
    set(handles.stba,'Enable','off')
    set(handles.checkMOV,'Enable','off')
    set(handles.checkOsc,'Enable','off')
    set(handles.checkReg,'Enable','off')
    set(handles.text14,'Enable','off')
    set(handles.text15,'Enable','off')
    set(handles.LPFtxt,'Enable','off')
    set(handles.HPFtxt,'Enable','off')
    set(handles.text21,'Enable','off')
    set(handles.addFunc,'String','Add FiltRegrSW Files')
else
    set(handles.addStru,'Enable','on')
    set(handles.FazFilt,'Enable','on')
    set(handles.LPF,'Enable','on')
    set(handles.HPF,'Enable','on')
    set(handles.stba,'Enable','on')
    set(handles.checkMOV,'Enable','on')
    set(handles.checkOsc,'Enable','on')
    set(handles.checkReg,'Enable','on')
    set(handles.text14,'Enable','on')
    set(handles.text15,'Enable','on')
    set(handles.LPFtxt,'Enable','on')
    set(handles.HPFtxt,'Enable','on')
    set(handles.text21,'Enable','on')
    set(handles.addFunc,'String','Add Functional Files')
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

function OutputFil_Callback(hObject, eventdata, handles)
function edit9_Callback(hObject, eventdata, handles)
function edit10_Callback(hObject, eventdata, handles)
function SizeVoxFunc_Callback(hObject, eventdata, handles)
function SizeMatFunc_Callback(hObject, eventdata, handles)
function numDyna_Callback(hObject, eventdata, handles)
function LPF_Callback(hObject, eventdata, handles)
function checkMOV_Callback(hObject, eventdata, handles)
function checkOsc_Callback(hObject, eventdata, handles)
function HPF_Callback(hObject, eventdata, handles)
function stba_Callback(hObject, eventdata, handles)

function LPF_CreateFcn(hObject, eventdata, handles)
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
function wsIn_Callback(hObject, eventdata, handles)

function wsIn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ROIinp_Callback(hObject, eventdata, handles)

function ROIinp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TRInp_Callback(hObject, eventdata, handles)

function TRInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SSInp_Callback(hObject, eventdata, handles)

function SSInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
