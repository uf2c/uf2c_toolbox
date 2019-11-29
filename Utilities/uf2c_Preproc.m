function uf2c_Preproc(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX,NativeS)
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
%
    UDV = uf2c_defaults('Preproc_uf2c');
    
    if ~exist('NativeS','var')
        NativeS = 0;
    end

    if isequal(NativeS,0) % MNI Space preprocessing
        spm_jobman('initcfg')
        clear flags
        flags.quality = UDV.RQua;
        flags.fwhm  = UDV.RSmo;
        flags.sep = UDV.RSep;
        flags.rtm = 1;
        flags.wrap = [0 1 0];
        flags.interp = UDV.RInte;
        spm_realign(char(file1),flags);
        
        clear flags
        flags.mask = 1;
        flags.mean  = 1;
        flags.interp = UDV.RResI;
        flags.which = UDV.RIRes;
        flags.wrap = [0 1 0];
        flags.prefix = 'r';
        spm_reslice(char(file1),flags);

        if STC
            spm_slice_timing(char(file1),1:Funpre.dat.dim(3),...
                round(Funpre.dat.dim(3)/2),...
                [(TR-(TR/Funpre.dat.dim(3)))/(Funpre.dat.dim(3)-1) TR-(TR-(TR/Funpre.dat.dim(3)))], 'a');
        end
        
        clear job
        job.ref = {[pathFunc,dirname,filesep, 'mean' fileFuncX]};
        job.source = {file2};
        job.other = {''};
        job.eoptions.cost_fun = 'nmi';
        job.eoptions.sep = UDV.CorSep;
        job.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        job.eoptions.fwhm = [7 7];
        spm_run_coreg(job);
        
        clear job flags Channel1 Channel2 Channel3 Channel4 Channel5 Channel6
        job.channel.vols = {file2};
        job.channel.biasreg = 0.001;
        job.channel.biasfwhm = 60;
        job.channel.write = [0 1];
            Channel1.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,1']};
            Channel1.ngaus = 1;
            Channel1.native = UDV.SegSN;
            Channel1.warped = UDV.SegSW;
        job.tissue(1) = Channel1;
            Channel2.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,2']};
            Channel2.ngaus = 1;
            Channel2.native = UDV.SegSN;
            Channel2.warped = UDV.SegSW;
        job.tissue(2) = Channel2;
            Channel3.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,3']};
            Channel3.ngaus = 2;
            Channel3.native = UDV.SegSN;
            Channel3.warped = UDV.SegSW;
        job.tissue(3) = Channel3;
            Channel4.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,4']};
            Channel4.ngaus = 3;
            Channel4.native = [0 0];
            Channel4.warped = [0 0];
        job.tissue(4) = Channel4;
            Channel5.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,5']};
            Channel5.ngaus = 4;
            Channel5.native = [0 0];
            Channel5.warped = [0 0];
        job.tissue(5) = Channel5;
            Channel6.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,6']};
            Channel6.ngaus = 2;
            Channel6.native = [0 0];
            Channel6.warped = [0 0];
        job.tissue(6) = Channel6;
        job.warp.mrf = 1;
        job.warp.cleanup = 1;
        job.warp.reg = [0 0.001 0.5 0.05 0.2];
        job.warp.affreg = 'mni';
        job.warp.fwhm = 0;
        job.warp.samp = UDV.SegSep;
        job.warp.write = [0 1];
        spm_preproc_run(job);
        clear job Channel1 Channel2 Channel3 Channel4 Channel5 Channel6

        clear job
        job.subj.def = {[pathFunc,dirname,filesep,'y_',fileStruX]};
        job.subj.resample = {[pathFunc,dirname,filesep,'m',fileStruX]};
        if BounB
            job.woptions.bb = [-78 -112  -50
                                78   74   85];
        else
            job.woptions.bb = UDV.MNIBB;
        end
        job.woptions.vox = [Svs(1) Svs(2) Svs(3)];
        job.woptions.interp = 1;
        job.woptions.prefix = 'w';
        spm_run_norm(job);

        if STC
            clear file1
            for tt = 1:Funpre.dat.dim(4)
                file1{tt,1} = [pathFunc,dirname,filesep,'a',fileFuncX,',',sprintf('%d',tt)];
            end
        end
        
        clear job
        job.subj.def = {[pathFunc,dirname,filesep,'y_',fileStruX]};
        job.subj.resample = file1;
        if BounB
            job.woptions.bb = [-78 -112  -50
                                78   74   85];
            job.woptions.vox = [fvs(1) fvs(2) fvs(3)];
        else
            job.woptions.bb = UDV.MNIBB;
            job.woptions.vox = UDV.MNIFVS;
        end
        job.woptions.interp = 4;
        job.woptions.prefix = 'w';
        spm_run_norm(job);

        if ~STC
            movefile([pathFunc,dirname,filesep,'w',fileFuncX],...
                [pathFunc,dirname,filesep,'sw',fileFuncX],'f');
        else
            movefile([pathFunc,dirname,filesep,'wa',fileFuncX],...
                [pathFunc,dirname,filesep,'sw',fileFuncX],'f');
        end
        
    else   % For Native Space Preprocessing
        
        clear matlabbatch
        
        clear flags
        flags.quality = UDV.RQua;
        flags.fwhm  = UDV.RSmo;
        flags.sep = UDV.RSep;
        flags.rtm = 1;
        flags.wrap = [0 0 0];
        flags.interp = UDV.RInte;
        spm_realign(char(file1),flags);
        
        clear flags
        flags.mask = 1;
        flags.mean  = 1;
        flags.interp = UDV.RResI;
        flags.which = UDV.RIResNS;
        flags.wrap = [0 0 0];
        flags.prefix = 'r';
        spm_reslice(char(file1),flags);

        for h = 1:size(file1,1)
            file1{h} = [pathFunc,dirname,filesep,'r',fileFuncX,',',num2str(h)];
        end

        if STC
            spm_slice_timing(char(file1),1:Funpre.dat.dim(3),...
                round(Funpre.dat.dim(3)/2),...
                [(TR-(TR/Funpre.dat.dim(3)))/(Funpre.dat.dim(3)-1) TR-(TR-(TR/Funpre.dat.dim(3)))],'a');
            for h = 1:size(file1,1)
                file1{h} = [pathFunc,dirname,filesep,'ar',fileFuncX,',',num2str(h)];
            end
        end

        clear job
        job.ref = {[pathFunc,dirname,filesep,'mean',fileFuncX]};
        job.source = {file2};
        job.other = {''};
        job.eoptions.cost_fun = 'nmi';
        job.eoptions.sep = UDV.CorSep;
        job.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        job.eoptions.fwhm = [7 7];
        spm_run_coreg(job);
        
        clear job flags Channel1 Channel2 Channel3 Channel4 Channel5 Channel6
        job.channel.vols = {file2};
        job.channel.biasreg = 0.001;
        job.channel.biasfwhm = 60;
        job.channel.write = [0 1];
            Channel1.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,1']};
            Channel1.ngaus = 1;
            Channel1.native = UDV.SegSNns;
            Channel1.warped = UDV.SegSWns;
        job.tissue(1) = Channel1;
            Channel2.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,2']};
            Channel2.ngaus = 1;
            Channel2.native = UDV.SegSNns;
            Channel2.warped = UDV.SegSWns;
        job.tissue(2) = Channel2;
            Channel3.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,3']};
            Channel3.ngaus = 2;
            Channel3.native = UDV.SegSNns;
            Channel3.warped = UDV.SegSWns;
        job.tissue(3) = Channel3;
            Channel4.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,4']};
            Channel4.ngaus = 3;
            Channel4.native = [0 0];
            Channel4.warped = [0 0];
        job.tissue(4) = Channel4;
            Channel5.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,5']};
            Channel5.ngaus = 4;
            Channel5.native = [0 0];
            Channel5.warped = [0 0];
        job.tissue(5) = Channel5;
            Channel6.tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,6']};
            Channel6.ngaus = 2;
            Channel6.native = [0 0];
            Channel6.warped = [0 0];
        job.tissue(6) = Channel6;
        job.warp.mrf = 1;
        job.warp.cleanup = 1;
        job.warp.reg = [0 0.001 0.5 0.05 0.2];
        job.warp.affreg = 'mni';
        job.warp.fwhm = 0;
        job.warp.samp = UDV.SegSep;
        job.warp.write = [0 0];
        spm_preproc_run(job);
        clear job Channel1 Channel2 Channel3 Channel4 Channel5 Channel6

        if STC
            movefile([pathFunc,dirname,filesep,'ar',fileFuncX],...
                [pathFunc,dirname,filesep,'sr',fileFuncX],'f');
        else
           movefile([pathFunc,dirname,filesep,'r',fileFuncX],...
                [pathFunc,dirname,filesep,'sr',fileFuncX],'f');
        end
    end
end


