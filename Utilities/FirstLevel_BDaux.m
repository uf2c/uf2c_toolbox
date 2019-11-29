function FirstLevel_BDaux(fullpath_folder,EPI_pref,task_name,onsets,duration,...
    TR,RP,contrast,pthresh,mult_cor,clust,FDRcInfe,glass_b,sl_view,T1_pref,SlicesN,Orie)
% FIRSTLEVEL_BDaux MATLAB code
%
% UF²C - User Friendly Functional Connectivity
% Raphael Fernandes Casseb
% University of Campinas, 2016
%
% Copyright (c) 2016, Raphael Fernandes Casseb
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
%
% This script is part of RS_ROI and runs first level analysis of task
% images.
% Images of a two-block-design are expected (rest and task). It runs the 
% following modules of the stats batch: fMRI model specification, model 
% estimation and contrast manager. It also saves a Glass Brain image of
% the selected contrast.
% Only Task>Rest and Rest>Task constrasts are automatically generated.

current_folder = pwd;

%%% Brunno Virus (Thresholds strings)
pthreshSTR = num2str(pthresh);
pthreshSTR = strrep(pthreshSTR,'.','');
clustSTR = num2str(clust);

if isequal(mult_cor,'none')
    mult_corSTR = 'Uncorr';
else
    mult_corSTR = mult_cor;
end
%%% 

cd(fullpath_folder)
saving_dir = {fullfile(fullpath_folder,['Stat_' pthreshSTR mult_corSTR '_c' clustSTR filesep])};

mkdir(saving_dir{1})



% Storing func and struc files names:
% -----------------------------------
EPI = dir([EPI_pref '*.nii']);

for i = 1:size(EPI,1) % For each EPI
    EPIaux = nifti(EPI(i,1).name);
    for t=1:EPIaux.dat.dim(4)
        eval(sprintf('EPI_%d{%d,1} = ''%s,%s'';',i,t,fullfile(fullpath_folder,EPI(i,1).name),num2str(t)));
    end
    clear EPIaux
end


% Moviment parameter files
% (There must be the same amount of EPIs)
% ---------------------------------------
rp_file = dir('rp_*.txt');
for i = 1:size(EPI,1) % For each EPI
    if ~RP % If rp files are not going to be used:
        eval(sprintf('rp_file_%d = '''';',i));
    else
        eval(sprintf('rp_file_%d = ''%s'';',i,fullfile(fullpath_folder,...
            rp_file(i,1).name)));
    end
end


if sl_view
    T1name  = dir([T1_pref '*.nii']);
    file_T1 = ([T1name.folder filesep T1name.name]);
end

% Contrast vector
% ---------------
if RP
    con1 = [];
    for cc = 1:size(EPI,1)
        RegTest = load(eval(sprintf('rp_file_%d',cc)));
        SizeRegT = size(RegTest,2);
        conTmp = ([1,zeros(1,SizeRegT)]);
        con1 = [con1,conTmp];
    end
    con2 = -con1;
else
    con1 = zeros(1,3*size(EPI,1));
    for cc=1:size(EPI,1)
        con1(1,cc*3-2) = 1;
    end
    con2 = -con1;
end

spm('defaults','fmri')
spm_jobman('initcfg')
clear matlabbatch

matlabbatch{1}.spm.stats.fmri_spec.dir = saving_dir;
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
if size(EPI,1)>1
    for iii = 1:size(EPI,1)
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).scans = ...
            eval(sprintf('EPI_%d',iii));
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.name = task_name;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.onset = onsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.duration = duration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).cond.orth = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).multi_reg = ...
            {eval(sprintf('rp_file_%d',iii))};
        matlabbatch{1}.spm.stats.fmri_spec.sess(iii).hpf = 128;
    end
else
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = EPI_1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = task_name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = onsets;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = duration;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rp_file_1};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
end

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [task_name '>Rest'];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = con1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['Rest>' task_name];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = con2;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

save(fullfile(saving_dir{1},'1st_level_batch'),'matlabbatch')
[output_list, ~] = spm_jobman('run',matlabbatch);

clear matlabbatch
%%% Brunno Virus (save thresholded SPM and Tiff for results
matlabbatch{1}.spm.stats.results.spmmat(1) = {[saving_dir{1} filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
matlabbatch{1}.spm.stats.results.conspec.extent = clust;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'tif';
matlabbatch{1}.spm.stats.results.write.tspm.basename = [task_name '_Posit_' pthreshSTR mult_corSTR '_c' clustSTR];
[output_list, ~] = spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.results.spmmat(1) = {[saving_dir{1} filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
matlabbatch{1}.spm.stats.results.conspec.extent = clust;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = 'tif';
matlabbatch{1}.spm.stats.results.write.tspm.basename = [task_name '_Negat_' pthreshSTR mult_corSTR '_c' clustSTR];
spm_jobman('run',matlabbatch)
try
    [~] = uf2c_AnatDescrip(saving_dir{1},[mult_corSTR,'_AnatDescrip.txt'],...
        {[saving_dir{1},filesep,'spmT_0001_',task_name,'_Posit_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.nii'];...
        [saving_dir{1},filesep,'spmT_0002_',task_name,'_Negat_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.nii']},...
        'map');
catch
    fprintf('Anatomical Description tool: The functional file resolution does not match our anatomical atlas')
end

if FDRcInfe
    
    FDRcV_Pos = round(output_list{1, 1}.TabDatvar.ftr{5,2}(4));
    
    clear matlabbatch
    matlabbatch{1}.spm.stats.results.spmmat(1) = {[saving_dir{1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
    matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
    matlabbatch{1}.spm.stats.results.conspec.extent = FDRcV_Pos;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.print = 'tif';
    matlabbatch{1}.spm.stats.results.write.tspm.basename = [task_name '_Posit_Init' pthreshSTR mult_corSTR '_FDRclc_' num2str(FDRcV_Pos) 'vox'];
    spm_jobman('run',matlabbatch)
    
    clear matlabbatch
    matlabbatch{1}.spm.stats.results.spmmat(1) = {[saving_dir{1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts = 2;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
    matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
    matlabbatch{1}.spm.stats.results.conspec.extent = FDRcV_Pos;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.print = 'tif';
    matlabbatch{1}.spm.stats.results.write.tspm.basename = [task_name '_Negat_Init' pthreshSTR mult_corSTR '_FDRclc_' num2str(FDRcV_Pos) 'vox'];
    spm_jobman('run',matlabbatch)
    try    
        [~] = uf2c_AnatDescrip(saving_dir{1},'FDR_AnatDescrip.txt',...
            {[saving_dir{1},'spmT_0001_',task_name,'_Posit_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.nii'];...
            [saving_dir{1},'spmT_0002_',task_name,'_Negat_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.nii']},...
            'map');
    catch
        fprintf('Anatomical Description tool: The functional file resolution does not match our anatomical atlas')
    end
end

if glass_b
    % Generating glass brain image\
    xSPM.swd   = saving_dir{1};             % SPM.mat folder
    xSPM.u     = pthresh;                % Threshold
    xSPM.Im    = [];                     % Masking
    xSPM.Ex    = [];                     % exclusive or inclusive
    xSPM.k     = clust;                  % Extent (voxels)
    xSPM.units = {'mm','mm','mm'};       % Data type: Volumetric (2D/3D)
    xSPM.title     = 'Task>Rest';        % Results title
    xSPM.Ic        = contrast;           % Contrast(s)
    xSPM.thresDesc = mult_cor;           % Threshold type: 'none', 'FWE' (must be inside the for-loop)

    [SPM,xSPM] = spm_getSPM(xSPM);
    
    MIP = spm_mip(xSPM.Z,xSPM.XYZmm,xSPM.M);
    imwrite(MIP,gray(64),fullfile(saving_dir{1},'Glass_brain.png'),'png');
    
end

if sl_view
    namePos = [saving_dir{1},filesep,'spmT_0001_',task_name,'_Posit_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.nii'];

    nameNeg = [saving_dir{1},filesep,'spmT_0002_',task_name,'_Negat_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.nii'];

    SliceView_uf2c(file_T1,namePos,0.1,0.1,Orie,'hot',...
        SlicesN,'',['SliceView_spmT_0001_',task_name,'_Posit_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.png'],...
            saving_dir{1})
    SliceView_uf2c(file_T1,nameNeg,0.1,0.1,Orie,'winter',...
        SlicesN,'',['SliceView_spmT_0002_',task_name,'_Negat_',pthreshSTR,mult_corSTR,'_c',clustSTR,'.png'],...
            saving_dir{1})
    if FDRcInfe    
        namePos2 = [saving_dir{1},filesep,'spmT_0001_',task_name,'_Posit_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.nii'];

        nameNeg2 = [saving_dir{1},filesep,'spmT_0002_',task_name,'_Negat_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.nii'];

        SliceView_uf2c(file_T1,namePos2,0.1,0.1,Orie,'hot',...
            SlicesN,'',['SliceView_spmT_0001_',task_name,'_Posit_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.png'],...
                saving_dir{1})
        SliceView_uf2c(file_T1,nameNeg2,0.1,0.1,Orie,'winter',...
            SlicesN,'',['SliceView_spmT_0002_',task_name,'_Negat_Init',pthreshSTR,mult_corSTR,'_FDRclc_',num2str(FDRcV_Pos),'vox.png'],...
                saving_dir{1})
    end
end

    
clear saving_dir Regs_aux Mov_Reg EPI file_EPI scans MIP

cd(current_folder)

end