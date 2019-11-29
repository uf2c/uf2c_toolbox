function [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = mask_globals_uf2c(fvs,pathFunc,dirname,fileFunc,fileStru,SPMbb,FI,NativeS)
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

% This function mask the functional data and extract globals signals for
% regressions.
% the function assumes that the "dirname" folder contains normalized
% segmented mpas of the structural image (fileStru): wc1*fileStru, 
% wc2*fileStru, wc3*fileStru.
%
% The imputs are:
%
% fvs: vector 1x3 with the functional image voxel sizes / vector, double
% pathFunc: directory where the functional image folder is/ string
% dirname: name of the folder where the functional image is/ string
% fileFunc: name of the functional image (no path)/ string
% fileStru: name of the structural image (no path)/ string
% SPMbb: binary variable indicating the bounding box size (UF²C default: 0) / double
% FI: binary variable indicating if it is a Functional Interactivity
%   experiment: 0 or 1; This affects the brain masks thresholds
%
% Ex: [sizeFunc,finalEPI,MeanWM,MeanCSF] = mask_globals_uf2c([3,3,3],'C:\Users\fMRI\Desktop\Tests','Vol1','Vol1_fMRI.nii','Vol1_T1.nii',0,0)
%
% The outputs are:
%
% sizeFunc: The matrix sizes of the final 4D data.
% finalEPI:4D functional data double variable. Masked by the interpolated
% 	grey matter mask.
% MeanWM: vector (1x4th dim) with the average white matter signals time
% series
% MeanCSF: vector (1x4th dim) with the average CSF signals time
% series
% Func: struct with the functional data header
%

    if ~exist('NativeS','var')
        NativeS = 0;
    end

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Extracting fMRI global signals\n');
    
    if isequal(NativeS,0) % MNI Space preprocessing
        
        Func = nifti([pathFunc,dirname,filesep,'sw',fileFunc]);
        Struct = nifti([pathFunc,dirname,filesep,'wm',fileStru]);

        matrix4 = Func.dat(:,:,:,:);
        sizeFunc = size(matrix4);

        template = matrix4;
        template(:,:,:) = 1;
        templ = Func;
        templ.dat.dim = [sizeFunc(1) sizeFunc(2) sizeFunc(3)];
        templ.dat.fname = [pathFunc,dirname,filesep,'Template.nii'];
        templ.dat.dtype = 'FLOAT32-LE';
        templ.dat(:,:,:) = template;
        create(templ)
        clear template

        %%%%%%%%%%%%%%%%%%%%%%% Interp segmented images %%%%%%%%%%%%%%%%%%%%%%%%%%
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %GM image
                                            [pathFunc,dirname,filesep,'mwc1',fileStru]};
        matlabbatch{1}.spm.util.imcalc.output = 'c1interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.2';% Apply threshod to the image. Removes voxels with less than 25% of chance to be GM.
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc,dirname,filesep,'c1interp.nii']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
        
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch)

        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'sc1interp.nii']};
        matlabbatch{1}.spm.util.imcalc.output = 'c1interpF.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.35';% Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)

        c1Nii = nifti([pathFunc,dirname,filesep,'c1interpF.nii']);
        c1mat = c1Nii.dat(:,:,:);

        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %WM image
                                            [pathFunc,dirname,filesep,'mwc2',fileStru]};
        matlabbatch{1}.spm.util.imcalc.output = 'c2interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)

        c2Nii = nifti([pathFunc,dirname,filesep,'c2interp.nii']);
        c2mat = c2Nii.dat(:,:,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%  To include the WM on the analysis mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[pathFunc,dirname,filesep,'c2interp.nii']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch)
    
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'sc2interp.nii']};
        matlabbatch{1}.spm.util.imcalc.output = 'c2interpF.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.35';% Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch)
    
        c2Nii = nifti([pathFunc,dirname,filesep,'c2interpF.nii']);
        c2mat = c2Nii.dat(:,:,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[pathFunc,dirname,filesep,'Template.nii']  %CSF image
                                                [pathFunc,dirname,filesep,'mwc3',fileStru]};
        matlabbatch{1}.spm.util.imcalc.output = 'c3interp.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {[pathFunc,dirname]};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2>0.4';
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

        for vx = 1:sizeFunc(4) 
%              finalEPI(:,:,:,vx) = c1mat.*matrix4(:,:,:,vx);% Apply Grey Matter mask to the original EPI to remove outside brain voxels
            finalEPI(:,:,:,vx) = ((c1mat+c2mat)>0).*matrix4(:,:,:,vx); % Alternative mask with GM and WM

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

        save([pathFunc,dirname,filesep,'MeanWM'],'MeanWM')
        save([pathFunc,dirname,filesep,'MeanCSF'],'MeanCSF')

        clear c1mat c2mat c3mat c1Nii c2Nii c3Nii

        fprintf('Done! ============== %s\r\n',spm('time'));
        
end