function varargout = GeneralPreProc(varargin)
%
% By default, the directories should to contains functional files and a structural
% image. You can add many folders. Folders with more than one functional
% file will be considered a mult-session expirement.
% Use always NIfTI format.
% KEEP IN THE FOLDER ONLY THE (RAW) NIfTI FILES OF INTEREST.
%
% Brunno Machado de Campos
% University of Campinas, 2017

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
                   'gui_OpeningFcn', @GeneralPreProc_OpeningFcn, ...
                   'gui_OutputFcn',  @GeneralPreProc_OutputFcn, ...
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


function GeneralPreProc_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = GeneralPreProc_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddFold_Callback(hObject, eventdata, handles)
global filesFold

filesFold = uipickfiles('Output','cell','Prompt','Add all subjects directories','REFilter','\');

if isequal(filesFold,0)
    set(handles.addtxt,'String','No folders added')
else
    ncases = size(filesFold,2);
    set(handles.addtxt,'String',[num2str(ncases) ' folder(s) added!'])
end

function runB_Callback(hObject, eventdata, handles)
global filesFold

set(handles.text9,'String','Running...')
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps for functionals
Realign = get(handles.realig,'Value');
SlcTimming = get(handles.SlcTim,'Value');
Correg = get(handles.coreg,'Value');
NormalFunc = get(handles.funcNorm,'Value');
SmoothFunc = get(handles.smooth,'Value');

% Steps for structurals
NormalStru = get(handles.checkbox7,'Value');

% Normal. Standard
MNIstd = get(handles.MNIDef,'Value');
SPM8std = get(handles.SPM8def,'Value');
SPM12std = get(handles.SPM12def,'Value');
if MNIstd
    NormSTD = 'MNI';
end
if SPM8std
    NormSTD = 'SPM8';
end
if SPM12std
    NormSTD = 'SPM12';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d Folder(s) added \r\n',size(filesFold,2))

for pn = 1:size(filesFold,2)
    
    fprintf('Starting forder %d of %d \n',pn,size(filesFold,2))
    
    cd(filesFold{pn})
    dirstru = dir(filesFold{pn});
    nfiles = (size(dirstru,1))-2;
    filesF = cell(nfiles,1);
    
    imgRR = getframe(GeneralPreProc);
    imwrite(imgRR.cdata, [filesFold{pn},filesep,'1-Your_Choices.png']);
    
    ExiFunc = 0;
    ExiStru = 0;

    for j = 3:nfiles+2
        if isequal(dirstru(j,1).isdir,1)
            filesF{j-2} = [];
        else
            filesF{j-2} = dirstru(j,1).name;
        end
    end
    
    filesF(cellfun(@isempty,filesF)) = [];
    filesF = sort(filesF);
    as = 1;
    sn = 0;
    for k = 1:size(filesF,1)
        [patch, fileN, exte] = fileparts(filesF{k});
        if isequal(exte,'.nii')
            teste = nifti(filesF{k});
            if isequal(size(teste.dat.dim,2),3)
                eval(sprintf('Stru = ''%s'';',filesF{k}))
                sn = sn+1;
                as = as-1;
            else
                eval(sprintf('EPI_%d = ''%s'';',as,filesF{k}))
            end
           as = as+1;
        end
    end
    
    if exist('EPI_1','var')
        fprintf('     %d functional images identified\n',as-1)
        ExiFunc = 1;
        for h = 1:(as-1)
            eval(sprintf('teste2 = nifti(EPI_%d);',h))
            for t=1:teste2.dat.dim(4)
                eval(sprintf('EPIf_%d{%d,1} = [EPI_%d '','' num2str(%d)];',h,t,h,t));
            end
        end

        a1 = 'EPIf_1';

        for ne = 2:(as-1)
            a1 = sprintf('%s,EPIf_%d',a1,ne);
        end

        FuncRes = nifti(EPI_1);
        FuncSls = FuncRes.dat.dim(3);
        FuncPix = [FuncRes.hdr.pixdim(2) FuncRes.hdr.pixdim(3) FuncRes.hdr.pixdim(4)];
        FuncPix =  double(FuncPix(1,:));

        TRvalue = FuncRes.timing.tspace;
        
        if isequal(TRvalue,0) || TRvalue<1
            TRvalue = inputdlg('Please, set the TR value in Seconds');
            TRvalue = str2num(TRvalue{1});
        end
        
    else
        fprintf('     No functional images identified\n')
        set(handles.text9,'String','Aborted...')
        drawnow
        return
    end
    
    if sn>1
        fprintf('     More than 1 structural (3D) images were identified\n')
        fprintf('     Only one is expected, aborting....\r\n')
        set(handles.text9,'String','Aborted...')
        drawnow
        return
    end
    
    if exist('Stru','var')
        fprintf('     Structural images identified\r\n')
        ExiStru = 1;
        StruRes = nifti(Stru);
        StruResPix = [StruRes.hdr.pixdim(2) StruRes.hdr.pixdim(3) StruRes.hdr.pixdim(4)];
        StruResPix =  double(StruResPix(1,:));
    else
        fprintf('     No structural image identified\r\n')
    end
    
    fprintf('     The Functional Image TR is %d sec\n',TRvalue)
    
    clear matlabbatch
    if ExiFunc
        if Realign
            
            if NormalFunc
                matlabbatch{1}.spm.spatial.realign.estwrite.data = eval(sprintf('{%s}',a1));
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                spm_jobman('run',matlabbatch)
                clear matlabbatch
            else
                matlabbatch{1}.spm.spatial.realign.estwrite.data = eval(sprintf('{%s}',a1));
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                spm_jobman('run',matlabbatch)
                clear matlabbatch
                
                for kk = 1:(as-1)
                    for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                        eval(sprintf('EPIf_%d{%d} = [''r'' EPIf_%d{%d}];',kk,fr,kk,fr))
                    end
                end
            end
        end
        
        if SlcTimming
            for kk = 1:(as-1)
                matlabbatch{1}.spm.temporal.st.scans = eval(sprintf('{EPIf_%d}',kk));
                matlabbatch{1}.spm.temporal.st.nslices = FuncSls;
                matlabbatch{1}.spm.temporal.st.tr = TRvalue;
                matlabbatch{1}.spm.temporal.st.ta = TRvalue-TRvalue/FuncSls;
                matlabbatch{1}.spm.temporal.st.so = 1:FuncSls;
                matlabbatch{1}.spm.temporal.st.refslice = ceil(FuncSls/2);
                matlabbatch{1}.spm.temporal.st.prefix = 'a';
                                        
                spm_jobman('run',matlabbatch)
                clear matlabbatch

                for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                    eval(sprintf('EPIf_%d{%d} = [''a'' EPIf_%d{%d}];',kk,fr,kk,fr))
                end

            end
        end
    end
    
    if ExiFunc && ExiStru && Realign && Correg
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {['mean' EPI_1]};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {Stru};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                        
        spm_jobman('run',matlabbatch)
        clear matlabbatch
    end
    
    if ExiStru && NormalStru
        spmP =  which('spm');
        spmP = spmP(1:end-5);
        
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {Stru};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmP filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
        
        matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {['m' Stru]};
        
        switch NormSTD
            case 'MNI'
                matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-90  -126   -72
                                                                           90    90   108];
                matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = StruResPix;
            case 'SPM8'
                matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -50
                                                                           78    76    85];
                matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = StruResPix;
            case 'SPM12'
                matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112  -70
                                                                           78   76   85];
                matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = StruResPix;
        end
        matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';
                
        spm_jobman('run',matlabbatch)
        clear matlabbatch
        
        Stru = ['wm' Stru];
        
        if ExiFunc && NormalFunc && Correg
            if SlcTimming
                for kk = 1:(as-1)
                    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {['y_' Stru(3:end)]};
                    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = eval(sprintf('EPIf_%d',kk));
                    
                    switch NormSTD
                        case 'MNI'
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90  -126   -72
                                                                                       90    90   108];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                        case 'SPM8'
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -50
                                                                                          78    76    85];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = FuncPix;
                        case 'SPM12'
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112  -70
                                                                                           78   76   85];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = FuncPix;
                    end
                    FinalVoxSize = matlabbatch{1}.spm.spatial.normalise.write.woptions.vox;
                    
                    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
                    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
                
                    spm_jobman('run',matlabbatch)
                    clear matlabbatch
                
                    for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                        eval(sprintf('EPIf_%d{%d} = [''w'' EPIf_%d{%d}];',kk,fr,kk,fr))
                    end
                    
                end
            else
                for kk = 1:(as-1)
                    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {['y_' Stru(3:end)]};
                    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = eval(sprintf('EPIf_%d',kk));

                    switch NormSTD
                        case 'MNI'
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90  -126   -72
                                                                                       90    90   108];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                        case 'SPM8'
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -50
                                                                                       78    76    85];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = FuncPix;
                        case 'SPM12'
%                             matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112  -70
%                                                                                        78   76   85];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -112  -90
                                                                                           90   90   90];
                            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = FuncPix;
                    end
                    FinalVoxSize = matlabbatch{1}.spm.spatial.normalise.write.woptions.vox;
                    
                    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
                    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
                
                    spm_jobman('run',matlabbatch)
                    clear matlabbatch
                
                    for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                        eval(sprintf('EPIf_%d{%d} = [''w'' EPIf_%d{%d}];',kk,fr,kk,fr))
                    end
                end
            end
        end
    else
        if ExiFunc && NormalFunc
            spmP =  which('spm');
            spmP = spmP(1:end-5);
            for kk = 1:(as-1)
                matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {eval(sprintf('EPIf_%d{1}',kk))};
                matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = eval(sprintf('EPIf_%d',kk));
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {[spmP '\tpm\TPM.nii']};
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
                matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
             
                switch NormSTD
                    case 'MNI'
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-90  -126   -72
                                                                                      90    90   108];
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
                    case 'SPM8'
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78  -112   -50
                                                                                      78    76    85];
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = FuncPix;
                    case 'SPM12'
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112  -70
                                                                                      78   76   85];
                        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = FuncPix;
                end
                
                FinalVoxSize = matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox;
                
                matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
                matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
                
                spm_jobman('run',matlabbatch)
                clear matlabbatch
                
                for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                    eval(sprintf('EPIf_%d{%d} = [''w'' EPIf_%d{%d}];',kk,fr,kk,fr))
                end
            end
        end
    end
    
    if ExiFunc && SmoothFunc
        for kk = 1:(as-1)
            matlabbatch{1}.spm.spatial.smooth.data = eval(sprintf('EPIf_%d',kk));
            if FuncPix(3)>=3
                matlabbatch{1}.spm.spatial.smooth.fwhm = 2.*FuncPix;
            else
                matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
            end
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            
            spm_jobman('run',matlabbatch)
            clear matlabbatch
            for fr = 1:size(eval(sprintf('EPIf_%d',kk)),1)
                eval(sprintf('EPIf_%d{%d} = [''s'' EPIf_%d{%d}];',kk,fr,kk,fr))
            end

        end
    end
    ComplexRegressors = 1;
    
    if ComplexRegressors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot and save the motion series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mkdir([filesFold{pn},filesep,'Motion_Controls_Params'])
        clear RPs
        for rpsI = 1:(as-1)
            tmptrp = eval(sprintf('EPI_%d',rpsI));
            tmptrp = [tmptrp(1:end-3) 'txt'];
            RPs{rpsI} = sprintf('%s%srp_%s',filesFold{pn},filesep,tmptrp);
        end

        [Mot_fig,HDvV,HTvV] = uf2c_plot_motion(RPs,'on');
        imgRR = getframe(Mot_fig);
        imwrite(imgRR.cdata,[filesFold{pn},filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot.png']);
        saveas(Mot_fig,[filesFold{pn},filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot'],'fig')
        close(Mot_fig)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScreSize = get(0,'screensize');
        ScreSize = ScreSize(3:end);
        
        for rpsI = 1:(as-1)
            tmptEPIv = eval(sprintf('EPI_%d',rpsI));
            rpVez = [filesFold{pn},filesep,'rp_',tmptEPIv(1:end-3),'txt'];
            tmptEPIv = eval(sprintf('EPIf_%d',rpsI));
            swaVez = [filesFold{pn},filesep,tmptEPIv{1}(1:end-2)];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BRAMILA's Framewise Displacemente (FD)
            % Code from: Brain and Mind Lab at Aalto University
            % Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also 
            % Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Quantifying Framewise Displacemente (FD)\r\n');

            cfg.motionparam = rpVez;
            cfg.prepro_suite = 'spm';
            cfg.radius = 50;
            [FDts,~] = bramila_framewiseDisplacement(cfg);
            
            figuFD = figure('Visible','off');
            set(figuFD,'Name','Framewise Displacemente (FD) TS',...
                'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
                'Color',[1 0.94 0.86]);
            
            BinVetFDMEAN = double(FDts>0.75);
            BinVetFDMEAN = BinVetFDMEAN.*FDts;
            plot(FDts);
            hold on
            
            plot(ones(1,numel(FDts)).*0.75);
            area(BinVetFDMEAN,'FaceColor','r')
            ylabel('FD (mm)')
            xlabel('Time Points')
            title('Average Framewise Displacement TS','FontSize', 14);
            drawnow
            imgRR = getframe(figuFD);
            imwrite(imgRR.cdata, [filesFold{pn},filesep,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS_EPI' num2str(rpsI) '.png']);
            saveas(figuFD,[filesFold{pn},filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS_EPI' num2str(rpsI)],'fig')
            close(figuFD)
            fprintf('Done! ============== %s\r\n',spm('time'));
            
            FDCensur = 1;
            if FDCensur
                fprintf('\r\n')
                fprintf('UF²C =============== %s\n',spm('time'));
                fprintf('Creating FD Temporal Mask\r\n');

                FD_TM = FDts>0.75;
                fprintf('Done! ============== %s\r\n',spm('time'));
            else
                FD_TM = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extracting Globals and Masking
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MaskEPI = 1;
            try
                FI = 0;
                if MaskEPI
                    [~,FEPI,~,MeanCSF,FuncStruTT] = mask_globals_uf2c(FinalVoxSize,filesFold{pn},'',swaVez,Stru(3:end),0,0,1,0);
                    FuncStruTT.dat.fname = swaVez;
                    FuncStruTT.dat(:,:,:,:) = FEPI;
                    create(FuncStruTT)
                else
                    [~,~,~,MeanCSF,~] = mask_globals_uf2c(FinalVoxSize,filesFold{pn},'',swaVez,Stru(3:end),0,0,1,0);
                end
            catch
                set(handles.status,'String','Error...')
                warndlg(sprintf('An error occured during subject %d globals and mask extraction. Check your data.',yy), 'Process Aborted')
                fclose('all');
                return
            end
            try
                close(Mot_fig)
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BRAMILA DVARS - Computes Derivative VARiance
            % Code from: Brain and Mind Lab at Aalto University
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            extGMvar = 1;
            extWMvar = 1;
            extCSFvar = 0;
        %     DVARthresV = str2num(get(handles.DVARthres,'String'));
            DVARthresV = 7.5;
            tmpstru = nifti(swaVez);

            if extGMvar || extWMvar || extCSFvar
                [dvarsC1,dvarsC2,~,MeanDvars] = bramila_dvars_uf2c(filesFold{pn},'',tmpstru.dat(:,:,:,:),extGMvar,extWMvar,extCSFvar,DVARthresV);
            else
                MeanDvars = [];
            end
            dvarCensu = 1;
            if dvarCensu
                fprintf('\r\n')
                fprintf('UF²C =============== %s\n',spm('time'));
                fprintf('Creating DVARs Temporal Mask\r\n');

                dvars_TM = zeros(tmpstru.dat.dim(4),1);
                dvars_TM = MeanDvars>7.5;
                fprintf('Done! ============== %s\r\n',spm('time'));
            else
                dvars_TM = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            RP_Mov = 1;
            CSFOscilation = 1;
            
            [Frgval,Regstr,figu,thicksL] = regrprep_uf2c(rpsI,RP_Mov,'','',rpVez,CSFOscilation,0,ScreSize,[],MeanCSF,tmpstru.dat.dim(4),0);
            if ~isempty(figu)
                imgRR = getframe(figu);
                imwrite(imgRR.cdata, [filesFold{pn},filesep,'Regression Matrix.tif']);
            end
            
            TM_DVARsV = 1;
            TM_FDV = 1;
            DVAR_tsV = 0;
            
            if TM_DVARsV || TM_FDV || DVAR_tsV
               [Frgval,figu,Censu] = Add_FD_DVAR_uf2c(figu,thicksL,Frgval,filesFold{pn},'',TM_DVARsV,TM_FDV,DVAR_tsV,MeanDvars,FD_TM,dvars_TM,tmpstru.dat.dim(4));
                if size(Frgval,2)>1
                    Regstr = [Regstr,'_tmpMask'];
                end
                if TM_DVARsV || TM_FDV
                    nOfSensu = size(Censu,2);
                end
            end
            
            close(figu)
            copyfile(rpVez,[rpVez(1:end-4) '_ORIGINAL.txt']);
            Frgval = Frgval(:,1:end-1);
            fidta = fopen(rpVez,'wt');
            for ii = 1:size(Frgval,1)
                fprintf(fidta,'%g\t',Frgval(ii,:));
                fprintf(fidta,'\n');
            end
            fclose(fidta)
            movefile([filesFold{pn},filesep,'Regression Matrix.tif'],[rpVez(1:end-4) '.png'],'f')
        end
    end
    
    for cler = 1:(k-1)
        clear(sprintf('EPIf_%d',cler));
        clear(sprintf('EPI_%d',cler));
    end
    
    delete([filesFold{pn} filesep  Stru(3:end-4) '_seg8.mat'])
    delete([filesFold{pn} filesep 'c1interp.nii'])
    delete([filesFold{pn} filesep 'sc2interp.nii'])
    delete([filesFold{pn} filesep 'c2interpReg.nii'])
    delete([filesFold{pn} filesep 'c3interpReg.nii'])
    delete([filesFold{pn} filesep 'sc1interp.nii'])
    delete([filesFold{pn} filesep 'c2interpF.nii'])
    delete([filesFold{pn} filesep 'c1interpF.nii'])
    delete([filesFold{pn} filesep 'c2interp.nii'])
    delete([filesFold{pn} filesep 'c3interp.nii'])
    delete([filesFold{pn} filesep 'Template.nii'])
    delete([filesFold{pn} filesep 'Template.mat'])
    clear Stru a1 as fileN FuncPix StruResPix filesF FuncRes StruRes
    
end
set(handles.text9,'String','Done!')
drawnow
fprintf('All DONE!\r\n');

function MNIDef_Callback(hObject, eventdata, handles)
if get(handles.MNIDef,'Value')
    set(handles.SPM12def,'Value',0)
    set(handles.SPM8def,'Value',0)
else
    set(handles.MNIDef,'Value',1)
end
    

function SPM12def_Callback(hObject, eventdata, handles)
if get(handles.SPM12def,'Value')
    set(handles.MNIDef,'Value',0)
    set(handles.SPM8def,'Value',0)
else
    set(handles.MNIDef,'Value',1)
end


function SPM8def_Callback(hObject, eventdata, handles)
if get(handles.SPM8def,'Value')
    set(handles.SPM12def,'Value',0)
    set(handles.MNIDef,'Value',0)
else
    set(handles.MNIDef,'Value',1)
end

function checkbox7_Callback(hObject, eventdata, handles)
if get(handles.checkbox7,'Value')
    set(handles.MNIDef,'Enable','on')
    set(handles.SPM8def,'Enable','on')
    set(handles.SPM12def,'Enable','on')
else
    if isequal(get(handles.funcNorm,'Value'),0)
        set(handles.MNIDef,'Enable','off')
        set(handles.SPM8def,'Enable','off')
        set(handles.SPM12def,'Enable','off')
    end
end

function funcNorm_Callback(hObject, eventdata, handles)
if get(handles.funcNorm,'Value')
    set(handles.MNIDef,'Enable','on')
    set(handles.SPM8def,'Enable','on')
    set(handles.SPM12def,'Enable','on')
else
    if isequal(get(handles.checkbox7,'Value'),0)
        set(handles.MNIDef,'Enable','off')
        set(handles.SPM8def,'Enable','off')
        set(handles.SPM12def,'Enable','off')
    end
end

function pushbutton5_Callback(hObject, eventdata, handles)
texteADD = {'You can add multiple folders with a set of images inside. Each folder could contain functional (4D) imageS';...
            'and a structural (3D) image, just functional files, or yet, just  a structural image. You  can add folders';...
            'with multiples functional images, but always with only one structural.';...
            'The toolbox will considerate a folders with more than one functional (4D) files as a multsession experiment.';...
            'All *.nii files inside each folder will me processed, so, keep just files of interest inside then.'};
set(handles.text7,'String',texteADD)

function pushbutton3_Callback(hObject, eventdata, handles)
texteSTD = {'-- The "MNI" option will result in normalized images with a matrix of 91 109 91 and voxel sizes of 2x2x2mm³';...
            '-- The "SPM12" option will result in a matrix of 53 63 52 for normalized FUNCTIONAL (4D) images and ';...
            '157 189 156 for normalized STRUCTURAL (3D) images. The voxel sizes will be kept equal to the raw images';...
            '-- The "SPM8" option will result in a matrix of 53 63 46 for FUNCTIONAL(4D) normalized images and';...
            '157 189 136 for STRUCTURAL(3D) normalized images. The voxel sizes will be kept equal to the raw images';...
            '-- For custom parameters, feel free to edit the code.'};
set(handles.text7,'String',texteSTD)

function realig_Callback(hObject, eventdata, handles)

function SlcTim_Callback(hObject, eventdata, handles)

function coreg_Callback(hObject, eventdata, handles)

function smooth_Callback(hObject, eventdata, handles)

function raphabtm_Callback(hObject, eventdata, handles)
