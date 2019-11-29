function UDV = uf2c_defaults(Modality)
if exist('Modality','var')
    switch Modality
%% UF²C Defaults
% UF²C parameters defined that is not customizable with the user graphical
% interfaces. Modify these parameters with caution.

        case 'FuncCon' 
%% 1. UF²C Defaults or defined parameters for 'FuncConn' Function

% 1.1. Parpool = Number of available cores minus the number of  
% cores to reserve (NC2R) to other machine tasks (default: 1)
            UDV.NC2R = 1;

% 1.2. Motion figure visible
            UDV.MFVstr = 'on';
            
% 1.3. FD brain radius
            UDV.FDBR = 50;
            
% 1.4. Functional Interactivity Study
            UDV.FI = 0;

% 1.5. Include WI on final image
            UDV.WMI = 0;

% 1.6. Native Space
            UDV.NativeS = 0;

% 1.7. Smooth factor: default is 3* voxel sizes
            UDV.SmoFac = 3;

% 1.8. ROI mask Outlier rigor. Default: 'severe'
            UDV.RMOR = 'severe';

% 1.9. SliceView Output
% 1.9.1 Overlay Threshold (r-score): 0.2
            UDV.SVOT = 0.2;
% 1.9.2 Underlay Threshold (intensity): 0.1
            UDV.SVUT = 0.1;
% 1.9.3 Plan. Default 'axial'
            UDV.SVPS = 'axial';
% 1.9.4 Number of slices. Default 60
            UDV.SVNoS = 60;
            
% 1.10. Maximum size of the Subject Dir Name
            UDV.MSDN = 30;
            
        case 'TS_CrossCor' 
%% 2. UF²C Defaults or defined parameters for 'TS_CrossCor' Function

% 2.1. Parpool = Number of available cores minus the number of  
% cores to reserve (NC2R) to other machine tasks (default: 1)
            UDV.NC2R = 1;
            
% 2.2. Motion figure visible
            UDV.MFVstr = 'on';
            
% 2.3. FD brain radius
            UDV.FDBR = 50;
            
% 2.4. Functional Interactivity Study
            UDV.FI = 0;

% 2.5. Include WI on final image
            UDV.WMI = 0;

% 2.6. Native Space
            UDV.NativeS = 0;

% 2.7. Smooth factor: default is 3* voxel sizes
            UDV.SmoFac = 3;

% 2.8. ROI mask Outlier rigor. Default: 'severe'
            UDV.RMOR = 'severe';

% 2.9. Maximum size of the Subject Dir Name
            UDV.MSDN = 30;
            
% 2.10. Significancy Threshold for the Average Circular Connectime Map
            UDV.STCCM = 0.05;
            
        case 'Preproc_RMC'
%% 3. UF²C Defaults or defined parameters for 'Preproc_RMC' function    

% 3.1. Parpool = Number of available cores minus the number of  
% cores to reserve (NC2R) to other machine tasks (default: 1)
            UDV.NC2R = 1;
            
% 3.2. Motion figure visible
            UDV.MFVstr = 'off';
            
% 3.3. FD brain radius
            UDV.FDBR = 50;
            
% 3.4. Functional Interactivity Study
            UDV.FI = 0;

% 3.5. Include WI on final image
            UDV.WMI = 0;

% 3.6. Native Space
            UDV.NativeS = 0;

% 3.7. Smooth factor: default is 3* voxel sizes
            UDV.SmoFac = 3;

% 3.9. Maximum size of the Subject Dir Name
            UDV.MSDN = 30;
            
        case 'JustPreproc'
%% 4. UF²C Defaults or defined parameters for 'JustPreproc' function

% 4.1. Parpool = Number of available cores minus the number of  
% cores to reserve (NC2R) to other machine tasks (default: 1)
            UDV.NC2R = 1;
            
% 4.2. Motion figure visible
            UDV.MFVstr = 'on';
            
% 4.3. FD brain radius
            UDV.FDBR = 50;
            
% 4.4. Functional Interactivity Study
            UDV.FI = 0;

% 4.5. Include WI on final image
            UDV.WMI = 0;

% 4.6. Native Space
            UDV.NativeS = 0;

% 4.7. Smooth factor: default is 3* voxel sizes
            UDV.SmoFac = 3;

% 4.8. ROI mask Outlier rigor. Default: 'severe'
            UDV.RMOR = 'severe';

% 4.9. Maximum size of the Subject Dir Name
            UDV.MSDN = 30;
            
        case 'Preproc_uf2c'
%% 5. UF²C Defaults or defined parameters for 'Preproc_uf2c' function    

% 5.1. Realign
    % 5.1.1 Quality. Default = 1
            UDV.RQua = 1;
    % 5.1.2 Separation. Default = 3
            UDV.RSep = 3;
    % 5.1.3 Smoothing. Default = 6
            UDV.RSmo = 5;
    % 5.1.4 Interpolation. Default = 2
            UDV.RInte = 2;
    % 5.1.5 Reslice Images. Default = [0 1]
            UDV.RIRes = [0 1];
    % 5.1.7 Reslice Images. NATIVE SPACE ANALYSIS Default = [2 1]
            UDV.RIResNS = [2 1];
    % 5.1.7 Reslice Interpolation. Default = 4
            UDV.RResI = 4;   

% 5.2. Corregister
    % 5.2.1 Separation. Default = [4 2]
            UDV.CorSep = [4 2];

% 5.3. Segmentation
    % 5.3.1 Save Native tissues: no ([0 0])
            UDV.SegSN = [0 0];
    % 5.3.2 Save Warped Tissues: no ([1 0])
            UDV.SegSW = [1 0];
    % 5.3.3 Sepration
            UDV.SegSep = 3;
            
% 5.4. Segmentation Case Native Space Analysis
    % 5.4.1 Save Native tissues: no ([0 0])
            UDV.SegSNns = [1 0];
    % 5.4.2 Save Warped Tissues: no ([1 0])
            UDV.SegSWns = [0 0];

% 5.5. Normalize
    % 5.5.1 MNI Bounding BOX. Default = [-90  -126   -72
    %                                   90    90   108]
            UDV.MNIBB = [-90  -126   -72
                          90    90   108];
    % 5.5.2 MNI Final Voxel sizes
            UDV.MNIFVS = [2 2 2];

        case 'NativeS_Conne' 
%% 6. UF²C Defaults or defined parameters for 'NativeS_Conne' Function

% 6.1. Parpool = Number of available cores minus the number of  
% cores to reserve (NC2R) to other machine tasks (default: 1)
            UDV.NC2R = 1;
            
% 6.2. Motion figure visible
            UDV.MFVstr = 'on';
            
% 6.3. FD brain radius
            UDV.FDBR = 50;
            
% 6.4. Functional Interactivity Study
            UDV.FI = 0;

% 6.5. Include WI on final image
            UDV.WMI = 0;

% 6.6. Native Space
            UDV.NativeS = 1;

% 6.7. Smooth factor: default is 2* voxel sizes
            UDV.SmoFac = 2;

% 6.8. ROI mask Outlier rigor. Default: 'severe'
            UDV.RMOR = 'severe';

% 6.9. Maximum size of the Subject Dir Name
            UDV.MSDN = 30;
            
% 6.10. Significancy Threshold for the Average Circular Connectime Map
            UDV.STCCM = 0.05;

        case 'imgthres_uf2c'
%% 7. UF²C Defaults or defined parameters for 'imgthres_uf2c' function 

% 7.1. Cutoff intensity. This is the threshold (in %)
% relative to the most frequenty image intensity
% (histogram peak, excluding zero). Default = 0.33 (30% of  
% the histogram peak intensity)
            UDV.ITP = 0.33;

        case 'mask_globals_uf2c'
%% 8. UF²C Defaults or defined parameters for 'mask_globals_uf2c' function

% 8.1. Apply threshod to the image. Removes voxels with less than 
% a XX% chance to be GM, on FI studies. Default(String) = '0.1'
            UDV.GMT1FI = '0.1';

% 8.2. Apply threshod to the image. Removes voxels with less than 
% a XX% chance to be GM. Default(String) = '0.1'
            UDV.GMT1 = '0.15';

% 8.3. Interpolation method for GM map. Default: -3 (3rd Degree Sinc)
            UDV.GMIMe = -3;

% 8.4. Smooth factor: default is 2* voxel sizes
            UDV.GMSFa = 2;

% 8.5. Apply threshod to the SMOOTHED GM image. Removes voxels  
% with less than a XX chance to be GM, on FI studies. Default
% (String) = '0.15'
            UDV.GMT2FI = '0.15';

% 8.6. Apply threshod to the SMOOTHED image. Removes voxels with 
% less than a XX% chance to be GM. Default(String) = '0.35'
            UDV.GMT2 = '0.35';

% 8.7. Apply threshold to WM tissue map to create the WM Mask 
% FOR REGRESSION. Removes voxels with less than a XX% chance 
% to be WM. Default: '0.95'.
            UDV.WMTR = '0.95';

% 8.8. Interpolation method for WM map. Default: -3 (3rd Degree Sinc)
            UDV.WMIMe = -3;

% 8.9. Apply threshold to WM tissue map to create the WM Mask
% FOR FMRI MASKING (to remove WM voxels). Default: '0.6'.
            UDV.WMTRm = '0.6';

% 8.10. Smooth factor for WM if you choose to include WM voxels 
% on final image: default is 2* voxel sizes
            UDV.WMSFa = 2;

% 8.11. Apply threshod to the SMOOTHED WM image if you choose  
% to include WM voxels on final image: Default(String) = '0.4'
            UDV.WMT2 = '0.4';

% 8.12. Apply threshold to CSF tissue map to create the CSF 
% Mask FORREGRESSION. Removes voxels with less than a XX% 
% chance to be CSF. Default: '0.95'.
            UDV.CSFTR = '0.95';

% 8.13. Interpolation method for CSF map. Default: -3 (3rd Degree Sinc)
            UDV.CSFIMe = -3;

% 8.14. Apply threshold to CSF tissue map to create the CSF 
% Mask FOR FMRI MASKING (to remove CSF voxels). Default: '0.6'.
            UDV.CSFTRm = '0.9';

% 8.15. Number of principal components extracted from WM 
% and CSF signals to include on the regression. Default: 5
% (first three) 
            UDV.NPCAs = 5;

        case 'bramila_dvars_uf2c'
%% 9. UF²C Defaults or defined parameters for 'bramila_dvars_uf2c' function

% 9.1. Intensity Factor. Default 1000
            UDV.FITF = 1000;

        case 'regrprep_uf2c'
%% 10. UF²C Defaults or defined parameters for 'regrprep_uf2c' function

        case 'Add_FD_DVAR_uf2c'
%% 11. UF²C Defaults or defined parameters for 'Add_FD_DVAR_uf2c' function  

        case 'filtregr_uf2c'
%% 12. UF²C Defaults or defined parameters for 'filtregr_uf2c' function     


%% End




    end
end
end
