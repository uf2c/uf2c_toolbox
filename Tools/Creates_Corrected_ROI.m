
FiltR = nifti('C:\Users\EEG-fMRI\Desktop\teste\BENEDITO_ROBERTO_DA_SILVA_EPI_SENSE_3_1\FiltRegrSW_BENEDITO_ROBERTO_DA_SILVA_EPI_SENSE_3_1.nii');
finalEPI = FiltR.dat(:,:,:,:);

ROI = nifti('C:\Users\EEG-fMRI\Desktop\teste\05-dorsalDMN-8.nii');
roi1 = ROI.dat(:,:,:);

Dx = FiltR.dat.dim(1);
Dy = FiltR.dat.dim(2);
Dz = FiltR.dat.dim(3);
Dt = FiltR.dat.dim(4);

% for fg = 1:size(roi1,3)
%     roi1(:,:,fg) = flipud(roi1(:,:,fg));
% end

for vx = 1:Dt; % get the ROI voxels time series. This step exclude from the ROI eventuals with matter voxels.
    final3D(:,:,:,vx) = roi1.*finalEPI(:,:,:,vx);
end

meanVet = squeeze(sum(sum(sum(final3D))))./nnz(mean(final3D,4));
FinalR = zeros(Dx,Dy,Dz);

for g = 1:Dx
    for j = 1:Dy
        for k = 1:Dz
            if sum(squeeze(final3D(g,j,k,:)))~= 0
                VetV = squeeze(final3D(g,j,k,:));
                xxX = corr(VetV,meanVet);
                FinalR(g,j,k) = xxX;
            end
        end
    end
end

FinalR(isnan(FinalR)) = 0;

FinalR_Z = 0.5.*(log(1+FinalR) - log(1-FinalR)); %Z transf

meanxxX = sum(sum(sum(FinalR_Z)))./nnz(FinalR_Z);

stdTmp1 = FinalR_Z-meanxxX;

stdTmp1(stdTmp1==-meanxxX)=0;
nnz(stdTmp1)

stdTmp1 = stdTmp1.^2;

stdxxX = sqrt(sum(sum(sum(stdTmp1)))./nnz(stdTmp1));

lowerCut = meanxxX-stdxxX;

FinalR_Z_Bin = FinalR_Z>=lowerCut;

SruF = ROI;
SruF.dat.fname = 'FinalROI.nii';
SruF.dat.dtype = 'FLOAT32-LE';
SruF.dat(:,:,:) = FinalR_Z;
create(SruF)

SruF = ROI;
SruF.dat.fname = 'FinalROI_Bin.nii';
SruF.dat.dtype = 'FLOAT32-LE';
SruF.dat(:,:,:) = FinalR_Z_Bin;
create(SruF)

