function [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = mask_globals_uf2c(fvs,pathFunc,dirname,fileFunc,fileStru,SPMbb,FI,WMI,NativeS)
%% UF²C - User Friendly Functional Connectivity
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
%%
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
%%

    UDV = uf2c_defaults('mask_globals_uf2c');

    if ~exist('WMI','var')
        WMI = 0;
    end

    if ~exist('NativeS','var')
        NativeS = 0;
    end

    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Extracting fMRI global signals\n');
    
    if isequal(NativeS,0) % MNI Space preprocessing
        
        try
            Func = nifti([pathFunc,dirname,filesep,'sw',fileFunc]);
        catch
            Func = nifti(fileFunc);
        end

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
        clear flags
        flags.dmtx = 0;
        flags.mask = 0;
        flags.interp = UDV.GMIMe;
        flags.dtype = 4;
        if FI
            equa = ['i1.*i2>',UDV.GMT1FI];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        else
            equa = ['i1.*i2>',UDV.GMT1];% Apply threshod to the image. Removes voxels with less than 35% of chance to be GM.
        end    
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %GM image
                    [pathFunc,dirname,filesep,'wc1',fileStru]},...
            [pathFunc,dirname,filesep,'c1interp.nii'],...
            equa,flags);

        % Smoothing
        clear job
        job.data = {[pathFunc,dirname,filesep,'c1interp.nii']};
        if SPMbb
            job.fwhm = [UDV.GMSFa*fvs(1) UDV.GMSFa*fvs(2) UDV.GMSFa*fvs(3)];
        else
            job.fwhm = [6 6 6];
        end
        job.dtype = 0;
        job.im = 0;
        job.prefix = 's';
        spm_run_smooth(job);

        clear flags
        flags.dmtx = 0;
        flags.mask = 0;
        flags.interp = UDV.GMIMe;
        flags.dtype = 4;
        if FI
            equa = ['i1>',UDV.GMT2FI];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        else
            equa = ['i1>',UDV.GMT2];% Apply threshod to the image. Removes voxels with less than 35% of chance to be GM.
        end    
        spm_imcalc({[pathFunc,dirname,filesep,'sc1interp.nii']},...
            [pathFunc,dirname,filesep,'c1interpF.nii'],...
            equa,flags);

        c1Nii = nifti([pathFunc,dirname,filesep,'c1interpF.nii']);
        c1mat = c1Nii.dat(:,:,:);

        % WM Mask for Regression
        flags.interp = UDV.WMIMe;
        equa = ['i1.*i2>',UDV.WMTR];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %WM image
                    [pathFunc,dirname,filesep,'wc2',fileStru]},...
            [pathFunc,dirname,filesep,'c2interpReg.nii'],equa,flags);

        % WM Mask for Regression
        c2NiiReg = nifti([pathFunc,dirname,filesep,'c2interpReg.nii']);
        c2matReg = c2NiiReg.dat(:,:,:);
        
        % WM Mask for MASKING
        equa = ['i1.*i2>',UDV.WMTRm];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %WM image
                    [pathFunc,dirname,filesep,'wc2',fileStru]},...
            [pathFunc,dirname,filesep,'c2interp.nii'],equa,flags);
        
        c2Nii = nifti([pathFunc,dirname,filesep,'c2interp.nii']);
        c2mat = c2Nii.dat(:,:,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%  To include the WM on the analysis mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if WMI
            clear job
            job.data = {[pathFunc,dirname,filesep,'c2interp.nii']};
            if SPMbb
                job.fwhm = [UDV.WMSFa*fvs(1) UDV.WMSFa*fvs(2) UDV.WMSFa*fvs(3)];
            else
                job.fwhm = [6 6 6];
            end
            job.dtype = 0;
            job.im = 0;
            job.prefix = 's';
            spm_run_smooth(job);
            
            clear flags
            flags.dmtx = 0;
            flags.mask = 0;
            flags.interp = UDV.WMIMe;
            flags.dtype = 4;
            equa = ['i1>',UDV.WMT2];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
            
            spm_imcalc({[pathFunc,dirname,filesep,'sc2interp.nii']},...
                [pathFunc,dirname,filesep,'c2interpF.nii'],...
                equa,flags);

            c2Nii = nifti([pathFunc,dirname,filesep,'c2interpF.nii']);
            c2mat = c2Nii.dat(:,:,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % CSF Mask for Regression
            clear flags
            flags.dmtx = 0;
            flags.mask = 0;
            flags.interp = UDV.CSFIMe;
            flags.dtype = 4;
            equa = ['i1.*i2>',UDV.CSFTR];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
            
            spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %CSF image
                        [pathFunc,dirname,filesep,'wc3',fileStru]},...
                [pathFunc,dirname,filesep,'c3interpReg.nii'],...
                equa,flags);

        c3NiiReg = nifti([pathFunc,dirname,filesep,'c3interpReg.nii']);
        c3matReg = c3NiiReg.dat(:,:,:);
        
        % CSF Mask for Masking
        equa = ['i1.*i2>',UDV.CSFTRm];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %CSF image
                    [pathFunc,dirname,filesep,'wc3',fileStru]},...
            [pathFunc,dirname,filesep,'c3interp.nii'],...
            equa,flags);
        
        c3Nii = nifti([pathFunc,dirname,filesep,'c3interp.nii']);
        c3mat = c3Nii.dat(:,:,:);

        finalEPI = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c2MapReg = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c3MapReg = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        
        if WMI
            MaskFinal = (c1mat+c2mat)-c3mat; % Alternative mask with GM and WM
            MaskFinal(MaskFinal<=0) = 0;
            MaskFinal(MaskFinal>0) = 1;
        else
            MaskFinal = (c1mat-c3mat)>0;% Apply Grey Matter mask to the original EPI to remove outside brain voxels
            MaskFinal(MaskFinal<=0) = 0;
            MaskFinal(MaskFinal>0) = 1;
        end
        for vx = 1:sizeFunc(4)
            finalEPI(:,:,:,vx) = MaskFinal.*matrix4(:,:,:,vx); % Alternative mask with GM and WM
            c2MapReg(:,:,:,vx) = c2matReg.*matrix4(:,:,:,vx); % Apply White Matter mask to the original EPI. Regress to WM variability
            c3MapReg(:,:,:,vx) = c3matReg.*matrix4(:,:,:,vx); % Apply CSF mask to the original EPI. regress to CSF variability
        end
        
        fprintf('Calculating WM PCAs\r\n');
        WMpcaRegs = uf2c_PCA(c2MapReg,UDV.NPCAs); %3
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        fprintf('Calculating CSF PCAs\r\n');
        CSFpcaRegs = uf2c_PCA(c3MapReg,UDV.NPCAs); %3
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        MeanWM = WMpcaRegs;
        MeanCSF = CSFpcaRegs;
        save([pathFunc,dirname,filesep,'PCA_WM'],'WMpcaRegs')
        save([pathFunc,dirname,filesep,'PCA_CSF'],'CSFpcaRegs')

        MaskFinalStr = c3Nii;
        MaskFinalStr.dat.fname = [pathFunc,dirname,filesep,'FinalMask.nii'];
        MaskFinalStr.dat(:,:,:) = MaskFinal;
        create(MaskFinalStr)
        
        clear c1mat c2mat c3mat c1Nii c2Nii c3Nii
        
    else  % Native Space Option
        
        Func = nifti(fileFunc);

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
        clear flags
        flags.dmtx = 0;
        flags.mask = 0;
        flags.interp = UDV.GMIMe;
        flags.dtype = 4;
        if FI
            equa = ['i1.*i2>',UDV.GMT1FI];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        else
            equa = ['i1.*i2>',UDV.GMT1];% Apply threshod to the image. Removes voxels with less than 35% of chance to be GM.
        end    

        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %GM image
                    [pathFunc,dirname,filesep,'c1',fileStru]},...
            [pathFunc,dirname,filesep,'c1interp.nii'],...
            equa,flags);

            clear job
            job.data = {[pathFunc,dirname,filesep,'c1interp.nii']};
            if SPMbb
                job.fwhm = [UDV.GMSFa*fvs(1) UDV.GMSFa*fvs(2) UDV.GMSFa*fvs(3)];
            else
                job.fwhm = [6 6 6];
            end
            job.dtype = 0;
            job.im = 0;
            job.prefix = 's';
            spm_run_smooth(job);
            
        if FI
            equa = ['i1>',UDV.GMT2FI];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        else
            equa = ['i1>',UDV.GMT2];% Apply threshod to the image. Removes voxels with less than 35% of chance to be GM.
        end    
        spm_imcalc({[pathFunc,dirname,filesep,'sc1interp.nii']},...
            [pathFunc,dirname,filesep,'c1interpF.nii'],...
            equa,flags);

        c1Nii = nifti([pathFunc,dirname,filesep,'c1interpF.nii']);
        c1mat = c1Nii.dat(:,:,:);

        % WM Mask for Regression
        flags.interp = UDV.WMIMe;
        equa = ['i1.*i2>',UDV.WMTR];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %WM image
                    [pathFunc,dirname,filesep,'c2',fileStru]},...
            [pathFunc,dirname,filesep,'c2interpReg.nii'],...
            equa,flags);
        
        c2NiiReg = nifti([pathFunc,dirname,filesep,'c2interpReg.nii']);
        c2matReg = c2NiiReg.dat(:,:,:);
        
        % WM Mask for MASKING
        equa = ['i1.*i2>',UDV.WMTRm];% Apply threshod to the image. Removes voxels with less than 15% chance to be GM.
        spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %WM image
                    [pathFunc,dirname,filesep,'c2',fileStru]},...
            [pathFunc,dirname,filesep,'c2interp.nii'],...
            equa,flags);
        
        c2Nii = nifti([pathFunc,dirname,filesep,'c2interp.nii']);
        c2mat = c2Nii.dat(:,:,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%%%%%  To include the WM on the analysis mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if WMI
            clear job
            job.data = {[pathFunc,dirname,filesep,'c2interp.nii']};
            if SPMbb
                job.fwhm = [UDV.WMSFa*fvs(1) UDV.WMSFa*fvs(2) UDV.WMSFa*fvs(3)];
            else
                job.fwhm = [6 6 6];
            end
            job.dtype = 0;
            job.im = 0;
            job.prefix = 's';
            spm_run_smooth(job);

            clear flags
            flags.dmtx = 0;
            flags.mask = 0;
            flags.interp = UDV.WMIMe;
            flags.dtype = 4;
            equa = ['i1>',UDV.WMT2]; %0.4 Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
            spm_imcalc({[pathFunc,dirname,filesep,'sc2interp.nii']},...
                [pathFunc,dirname,filesep,'c2interpF.nii'],...
                equa,flags);
            
            c2Nii = nifti([pathFunc,dirname,filesep,'c2interpF.nii']);
            c2mat = c2Nii.dat(:,:,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % CSF Mask for Regression
            clear flags
            flags.dmtx = 0;
            flags.mask = 0;
            flags.interp = UDV.CSFIMe;
            flags.dtype = 4;
            equa = ['i1.*i2>',UDV.CSFTR]; %0.4 Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
            spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  % CSF image
                       [pathFunc,dirname,filesep,'c3',fileStru]},...
                [pathFunc,dirname,filesep,'c3interpReg.nii'],...
                equa,flags);
        
        c3NiiReg = nifti([pathFunc,dirname,filesep,'c3interpReg.nii']);
        c3matReg = c3NiiReg.dat(:,:,:);
        
        % CSF Mask for Masking
            equa = ['i1.*i2>',UDV.CSFTRm]; %0.4 Apply threshod to the image. Removes voxels with more than 40% chance to not be GM.
            spm_imcalc({[pathFunc,dirname,filesep,'Template.nii']  %CSF image
                        [pathFunc,dirname,filesep,'c3',fileStru]},...
                [pathFunc,dirname,filesep,'c3interp.nii'],...
                equa,flags);
        
        c3Nii = nifti([pathFunc,dirname,filesep,'c3interp.nii']);
        c3mat = c3Nii.dat(:,:,:);

        finalEPI = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c2MapReg = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        c3MapReg = zeros(sizeFunc(1),sizeFunc(2),sizeFunc(3),sizeFunc(4)); % Prelocation
        
        if WMI
            MaskFinal = (c1mat+c2mat)-c3mat; % Alternative mask with GM and WM
            MaskFinal(MaskFinal<=0) = 0;
            MaskFinal(MaskFinal>0) = 1;
        else
            MaskFinal = (c1mat-c3mat)>0;% Apply Grey Matter mask to the original EPI to remove outside brain voxels
            MaskFinal(MaskFinal<=0) = 0;
            MaskFinal(MaskFinal>0) = 1;
        end
        
        for vx = 1:sizeFunc(4)
            
            finalEPI(:,:,:,vx) = MaskFinal.*matrix4(:,:,:,vx);
            
            c2MapReg(:,:,:,vx) = c2matReg.*matrix4(:,:,:,vx); % Apply White Matter mask to the original EPI. Regress to WM variability

            c3MapReg(:,:,:,vx) = c3matReg.*matrix4(:,:,:,vx); % Apply CSF mask to the original EPI. regress to CSF variability
        end
        
        fprintf('Calculating WM PCAs\r\n');
        WMpcaRegs = uf2c_PCA(c2MapReg,UDV.NPCAs); %3
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        fprintf('Calculating CSF PCAs\r\n');
        CSFpcaRegs = uf2c_PCA(c3MapReg,UDV.NPCAs); %3
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        MeanWM = WMpcaRegs;
        MeanCSF = CSFpcaRegs;
        save([pathFunc,dirname,filesep,'PCA_WM'],'WMpcaRegs')
        save([pathFunc,dirname,filesep,'PCA_CSF'],'CSFpcaRegs')

        MaskFinalStr = c3Nii;
        MaskFinalStr.dat.fname = [pathFunc,dirname,filesep,'FinalMask.nii'];
        MaskFinalStr.dat(:,:,:) = MaskFinal;
        create(MaskFinalStr)

        clear c1mat c2mat c3mat c1Nii c2Nii c3Nii

        fprintf('Done! ============== %s\r\n',spm('time'));
    end
end