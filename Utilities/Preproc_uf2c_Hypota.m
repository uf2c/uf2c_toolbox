function Preproc_uf2c_Hypota(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX,NativeS)
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
%
    if ~exist('NativeS','var')
        NativeS = 0;
    end

    if isequal(NativeS,0) % MNI Space preprocessing
        spm_jobman('initcfg')
        clear matlabbatch
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {file1};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 6;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 1 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 1 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

        spm_jobman('run',matlabbatch)

        if STC
            clear matlabbatch
            matlabbatch{1}.spm.temporal.st.scans = {file1};
            matlabbatch{1}.spm.temporal.st.nslices = Funpre.dat.dim(3);
            matlabbatch{1}.spm.temporal.st.tr = TR;
            matlabbatch{1}.spm.temporal.st.ta = TR-(TR/Funpre.dat.dim(3));
            matlabbatch{1}.spm.temporal.st.so = 1:Funpre.dat.dim(3);
            matlabbatch{1}.spm.temporal.st.refslice = round(Funpre.dat.dim(3)/2);
            matlabbatch{1}.spm.temporal.st.prefix = 'a';
            spm_jobman('run',matlabbatch)
        end

        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[pathFunc,dirname,filesep, 'mean' fileFuncX]};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {file2};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {file2};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pathFunc,dirname,filesep,'y_',fileStruX]};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pathFunc,dirname,filesep,'m',fileStruX]};
        if BounB
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112  -70
                                                                      78   76   85];
        else
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90  -126   -72
                                                                       90    90   108];
        end
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [Svs(1) Svs(2) Svs(3)];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
        spm_jobman('run',matlabbatch)

        if STC
            clear file1
            for tt = 1:Funpre.dat.dim(4)
                file1{tt,1} = [pathFunc,dirname,filesep,'a',fileFuncX,',',sprintf('%d',tt)];
            end
        end

        clear matlabbatch
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pathFunc,dirname,filesep,'y_',fileStruX]};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = file1;

        if BounB
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112  -70
                                                                       78   76   85];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [fvs(1) fvs(2) fvs(3)];
        else
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90  -126   -72
                                                                      90    90   108];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [fvs(1) fvs(2) fvs(3)];
        end

        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        spm_jobman('run',matlabbatch)

        clear file1
        if STC
            for tt = 1:Funpre.dat.dim(4)
                file1{tt,1} = [pathFunc,dirname,filesep,'wa',fileFuncX,',',sprintf('%d',tt)];
            end
        else
            for tt = 1:Funpre.dat.dim(4)
                file1{tt,1} = [pathFunc,dirname,filesep,'w',fileFuncX,',',sprintf('%d',tt)];
            end
        end

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = file1;

        if BounB
            matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        else
            matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 4.5];
        end

        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch)

        if STC
            fileFFunc = nifti([pathFunc,dirname,filesep,'swa',fileFuncX]);
            fileFFunc2 = fileFFunc;
            fileFFunc2.dat.fname = [pathFunc,dirname,filesep,'sw',fileFuncX];
            fileFFunc2.dat(:,:,:,:) = fileFFunc.dat(:,:,:,:);
            create(fileFFunc2)
            delete([pathFunc,dirname,filesep,'swa',fileFuncX]);
        end

        save([pathFunc,dirname,filesep,'Batch.mat'],'matlabbatch')
        
        
    else   % For Native Space Preprocessing
        
        spm_jobman('initcfg')
        clear matlabbatch
        
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {file1};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;
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
        
        for h = 1:size(file1,1)
            file1{h} = [pathFunc,dirname,filesep,'r',fileFuncX,',',num2str(h)];
        end

        if STC
            clear matlabbatch
            matlabbatch{1}.spm.temporal.st.scans = {file1};
            matlabbatch{1}.spm.temporal.st.nslices = Funpre.dat.dim(3);
            matlabbatch{1}.spm.temporal.st.tr = TR;
            matlabbatch{1}.spm.temporal.st.ta = TR-(TR/Funpre.dat.dim(3));
            matlabbatch{1}.spm.temporal.st.so = 1:Funpre.dat.dim(3);
            matlabbatch{1}.spm.temporal.st.refslice = round(Funpre.dat.dim(3)/2);
            matlabbatch{1}.spm.temporal.st.prefix = 'a';
            spm_jobman('run',matlabbatch)
            
            for h = 1:size(file1,1)
                file1{h} = [pathFunc,dirname,filesep,'ar',fileFuncX,',',num2str(h)];
            end

        end

        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[pathFunc,dirname,filesep,'mean',fileFuncX]};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {file2};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {file2};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[SPMdir2 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = file1;

        if BounB
            matlabbatch{1}.spm.spatial.smooth.fwhm = [2*fvs(1) 2*fvs(2) 2*fvs(3)];
        else
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
        end

        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch)

        if STC
            fileFFunc = nifti([pathFunc,dirname,filesep,'sar',fileFuncX]);
            fileFFunc2 = fileFFunc;
            fileFFunc2.dat.fname = [pathFunc,dirname,filesep,'sr',fileFuncX];
            fileFFunc2.dat(:,:,:,:) = fileFFunc.dat(:,:,:,:);
            create(fileFFunc2)
            delete([pathFunc,dirname,filesep,'sar',fileFuncX]);
        end

        save([pathFunc,dirname,filesep,'Batch.mat'],'matlabbatch')
    end

end


