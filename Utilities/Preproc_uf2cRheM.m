function Preproc_uf2cRheM(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX,NativeS)
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

        DirUF2C = which('uf2c');
        
        try
            pathTissue = [DirUF2C(1:end-6) 'Utilities' filesep 'RheMatlasfiles' filesep 'D99' filesep]; 
        catch
            pathTissue = [DirUF2C(1:end-13) 'Utilities' filesep 'RheMatlasfiles' filesep 'D99' filesep]; 
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
        matlabbatch{1}.spm.tools.oldseg.data = {file2};
        matlabbatch{1}.spm.tools.oldseg.output.GM = [1 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.WM = [1 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [1 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.biascor = 1;
        matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
        matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
                                                    [pathTissue 'c1D99_template.nii']
                                                    [pathTissue 'c2D99_template.nii']
                                                    [pathTissue 'c3D99_template.nii']
                                                    };
        matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2
                                                      2
                                                      2
                                                      4];
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'subj';
        matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
        matlabbatch{1}.spm.tools.oldseg.opts.warpco = 25;
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
        matlabbatch{1}.spm.tools.oldseg.opts.samp = 2;
        matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.tools.oldnorm.write.subj.matname = {[pathFunc,dirname,filesep,fileStruX(1:end-4),'_seg_sn.mat']};
        matlabbatch{1}.spm.tools.oldnorm.write.subj.resample = {[pathFunc,dirname,filesep,'m',fileStruX]};
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = [-34.2500  -49.7500  -27.5000
                                                               34.2500   36.7500   33.5000];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = [.5 .5 .5];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = 4;
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = 'w';        
        spm_jobman('run',matlabbatch)

        if STC
            clear file1
            for tt = 1:Funpre.dat.dim(4)
                file1{tt,1} = [pathFunc,dirname,filesep,'a',fileFuncX,',',sprintf('%d',tt)];
            end
        end

        clear matlabbatch
        matlabbatch{1}.spm.tools.oldnorm.write.subj.matname = {[pathFunc,dirname,filesep,fileStruX(1:end-4),'_seg_sn.mat']};
        matlabbatch{1}.spm.tools.oldnorm.write.subj.resample = file1;
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = [-34.2500  -49.7500  -27.5000
                                                               34.2500   36.7500   33.5000];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = [1 1 1];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = 4;
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = 'w';        
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
        matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
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
        
end


