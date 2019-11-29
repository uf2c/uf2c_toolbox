function [dvarsC1,dvarsC2,dvarsC3,MeanDvars] = uf2c_bramila_dvars(pathFunc,dirname,finalEPI,extGMvar,extWMvar,extCSFvar,DVARthresV)

UDV = uf2c_defaults('bramila_dvars_uf2c');

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Computing Derivative VARiance (DVARS)\r\n');

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

medI = median(finalEPI,4);

for t = 1:size(finalEPI,4)
    finalEPI(:,:,:,t) = finalEPI(:,:,:,t)./medI;
end

finalEPI(isinf(finalEPI)) = 0;
finalEPI(isnan(finalEPI)) = 0;
finalEPI = finalEPI.*UDV.FITF;

dvarsC1 = [];
dvarsC2 = [];
dvarsC3 = [];

Thresh = DVARthresV;

if extGMvar
    BramC1 = nifti([pathFunc,dirname,filesep,'c1interp.nii']);
    BramC1map = BramC1.dat(:,:,:);

    cfgC1.vol = finalEPI;
    cfgC1.plot = 0;
    cfgC1.mask = BramC1map;

    [dvarsC1,imgC1] = bramila_dvars(cfgC1);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_EPI',num2str(Nlo),'.mat'],'file')
            Nlo = Nlo+1;
        else
            save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_EPI',num2str(Nlo),'.mat'],'dvarsC1')
            Ok1 = 1;
        end
    end
    
    BinVetdvarC1 = double(dvarsC1>Thresh);
    BinVetdvarC1 = BinVetdvarC1.*dvarsC1;
    
    figuC1 = figure('Visible','off');
    set(figuC1,'Name','DerivativeVariance GM',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.3 ScreSize(1)*.3]),...
        'Color',[1 0.94 0.86]);
    imagesc(imgC1',[-2 2]);
    colormap(gray)
    title('Derivative Variance GM Masked Plot','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC1);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_EPI',num2str(Nlo),'.png']);
            Ok1 = 1;
        end
    end
    
%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM.png']);
    close(figuC1)

    figuC12 = figure('Visible','off');
    set(figuC12,'Name','DerivativeVariance GM TS',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
        'Color',[1 0.94 0.86]);
    plot(dvarsC1);
    hold on
    plot(ones(1,numel(dvarsC1)).*Thresh);
    area(BinVetdvarC1,'FaceColor','r')
    ylabel('DEVAR (%)')
    xlabel('Time Points')
    title('Average Derivative Variance GM Masked','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC12);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_avgTS_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_avgTS_EPI',num2str(Nlo),'.png']);
            saveas(figuC12,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_avgTS_EPI',num2str(Nlo)],'fig')
            Ok1 = 1;
        end
    end
    
%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_avgTS.png']);
%     saveas(figuC12,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_GM_avgTS'],'fig')
    
    delete(figuC12)
end

if extWMvar
    BramC2 = nifti([pathFunc,dirname,filesep,'c2interp.nii']);
    BramC2map = BramC2.dat(:,:,:);

    cfgC2.vol = finalEPI;
    cfgC2.plot = 0;
    cfgC2.mask = BramC2map;

    [dvarsC2,imgC2] = bramila_dvars(cfgC2);
%     save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM.mat'],'dvarsC2')
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_EPI',num2str(Nlo),'.mat'],'file')
            Nlo = Nlo+1;
        else
            save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_EPI',num2str(Nlo),'.mat'],'dvarsC1')
            Ok1 = 1;
        end
    end

    BinVetdvarC2 = double(dvarsC2>Thresh);
    BinVetdvarC2 = BinVetdvarC2.*dvarsC2;

    figuC2 = figure('Visible','off');

    set(figuC2,'Name','DerivativeVariance WM',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.3 ScreSize(1)*.3]),...
        'Color',[1 0.94 0.86]);
    imagesc(imgC2',[-2 2]);
    colormap(gray)
    title('Derivative Variance WM Masked Plot','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC2);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_EPI',num2str(Nlo),'.png']);
            Ok1 = 1;
        end
    end
%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM.png']);
    
    close(figuC2)

    figuC22 = figure('Visible','off');
    set(figuC22,'Name','DerivativeVariance WM TS',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
        'Color',[1 0.94 0.86]);
    plot(dvarsC2);
    hold on
    plot(ones(1,numel(dvarsC2)).*Thresh);
    area(BinVetdvarC2,'FaceColor','r')
    ylabel('DEVAR (%)')
    xlabel('Time Points')
    title('Average Derivative Variance WM Masked','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC22);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_avgTS_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_avgTS_EPI',num2str(Nlo),'.png']);
            saveas(figuC22,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_avgTS_EPI',num2str(Nlo)],'fig')
            Ok1 = 1;
        end
    end

%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_avgTS.png']);
%     saveas(figuC22,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_WM_avgTS'],'fig')
    delete(figuC22)
end

if  extCSFvar
    BramC3 = nifti([pathFunc,dirname,filesep,'c3interp.nii']);
    BramC3map = BramC3.dat(:,:,:);

    cfgC3.vol = finalEPI;
    cfgC3.plot = 0;
    cfgC3.mask = BramC3map;

    [dvarsC3,imgC3] = bramila_dvars(cfgC3);
%     save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF.mat'],'dvarsC3')
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_EPI',num2str(Nlo),'.mat'],'file')
            Nlo = Nlo+1;
        else
            save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_EPI',num2str(Nlo),'.mat'],'dvarsC1')
            Ok1 = 1;
        end
    end

    BinVetdvarC3 = double(dvarsC3>Thresh);
    BinVetdvarC3 = BinVetdvarC3.*dvarsC3;

    figuC3 = figure('Visible','off');

    set(figuC3,'Name','DerivativeVariance CSF',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.3 ScreSize(1)*.3]),...
        'Color',[1 0.94 0.86]);
    imagesc(imgC3',[-2 2]);
    colormap(gray)        
    title('Derivative Variance CSF Masked Plot','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC3);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_EPI',num2str(Nlo),'.png']);
            Ok1 = 1;
        end
    end

%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF.png']);
    close(figuC3)

    figuC32 = figure('Visible','off');
    set(figuC32,'Name','DerivativeVariance CSF TS',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
        'Color',[1 0.94 0.86]);
    plot(dvarsC3);
    hold on
    plot(ones(1,numel(dvarsC3)).*Thresh);
    area(BinVetdvarC3,'FaceColor','r')
    ylabel('DEVAR (%)')
    xlabel('Time Points')
    title('Average Derivative Variance CSF Masked','FontSize', 14);
    drawnow        
    imgRR = getframe(figuC32);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_avgTS_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_avgTS_EPI',num2str(Nlo),'.png']);
            saveas(figuC32,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_avgTS_EPI',num2str(Nlo)],'fig')
            Ok1 = 1;
        end
    end
    
%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_avgTS.png']);
%     saveas(figuC32,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_CSF_avgTS'],'fig')
    delete(figuC32)
end

    MeanDvars = zeros(size(dvarsC1,1),1);
    MeanDvarN = 0;
    if extGMvar
        MeanDvars = MeanDvars + dvarsC1;
        MeanDvarN = MeanDvarN + 1;
    end
    if extWMvar
        MeanDvars = MeanDvars + dvarsC2;
        MeanDvarN = MeanDvarN + 1;
    end
    if extCSFvar
        MeanDvars = MeanDvars + dvarsC3;
        MeanDvarN = MeanDvarN + 1;
    end
    
    MeanDvars = MeanDvars./MeanDvarN;

    BinVetdvarMEAN = double(MeanDvars>Thresh);
    BinVetdvarMEAN = BinVetdvarMEAN.*MeanDvars;

    figuMEAN = figure('Visible','off');
    set(figuMEAN,'Name','DerivativeVariance CSF TS',...
        'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
        'Color',[1 0.94 0.86]);
    plot(MeanDvars);
    hold on
    plot(ones(1,numel(MeanDvars)).*Thresh);
    area(BinVetdvarMEAN,'FaceColor','r')
    ylabel('DEVAR (%)')
    xlabel('Time Points')
    title('Average Derivative Variance AVERAGE tissues','FontSize', 14);
    drawnow        
    imgRR = getframe(figuMEAN);
    
    Ok1 = 0;
    Nlo = 1;
    while Ok1 == 0
        if exist([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_AVG_avgTS_EPI',num2str(Nlo),'.png'],'file')
            Nlo = Nlo+1;
        else
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_AVG_avgTS_EPI',num2str(Nlo),'.png']);
            saveas(figuMEAN,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_AVG_avgTS_EPI',num2str(Nlo)],'fig')
            Ok1 = 1;
        end
    end

%     imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_AVG_avgTS.png']);
%     saveas(figuMEAN,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'DerivativeVariance_AVG_avgTS'],'fig')
    
    delete(figuMEAN)

fprintf('Done! ============== %s\n\r',spm('time'));
end
