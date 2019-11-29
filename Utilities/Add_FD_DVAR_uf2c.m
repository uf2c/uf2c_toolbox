function [Frgval,figu,Censu] = Add_FD_DVAR_uf2c(figu,thicksL,Frgval,pathFunc,dirname,TM_DVARsV,TM_FDV,DVAR_tsV,MeanDvars,FD_TM,dvars_TM,dyn4)

    convVet = [1]; %Convoltion function

    if TM_DVARsV && TM_FDV

        TempMask = FD_TM + dvars_TM;

        TempMask(TempMask>1) = 1;
        TempMask = conv(convVet,double(TempMask));% convolution with a Convoltion function
        TempMask = TempMask(1:dyn4,1);
        TempMask(TempMask>1) = 1;

        TempMaskIdxs = find(TempMask);
        Censu = zeros(dyn4,numel(TempMaskIdxs));
        for i = 1:numel(TempMaskIdxs)
            Censu(TempMaskIdxs(i),i) = TempMask(TempMaskIdxs(i));
        end

        if ~isempty(Frgval)
            Frgval_tmp = Frgval(:,1:end-1);
            Frgval_tmp = [Censu,Frgval_tmp,ones(dyn4,1)];
            Frgval = Frgval_tmp;
        else
            Frgval_tmp = [Censu,ones(dyn4,1)];
            Frgval = Frgval_tmp;
        end
        thicksLCensor(1:size(Censu,2),1) = {'Censor'};
        thicksL = [thicksLCensor;thicksL];
    else
        if TM_FDV
            TempMask = FD_TM;

            TempMask(TempMask>1) = 1;
            TempMask = conv(convVet,double(TempMask));% convolution with a Convoltion function
            TempMask = TempMask(1:dyn4,1);
            TempMask(TempMask>1) = 1;

            TempMaskIdxs = find(TempMask);
            Censu = zeros(dyn4,numel(TempMaskIdxs));
            for i = 1:numel(TempMaskIdxs)
                Censu(TempMaskIdxs(i),i) = TempMask(TempMaskIdxs(i));
            end

            if ~isempty(Frgval)
                Frgval_tmp = Frgval(:,1:end-1);
                Frgval_tmp = [Frgval_tmp,Censu,ones(dyn4,1)];
                Frgval = Frgval_tmp;
            else
                Frgval_tmp = [Censu,ones(dyn4,1)];
                Frgval = Frgval_tmp;
            end
            thicksLCensor(1:size(Censu,2),1) = {'Censor'};
            thicksL = [thicksLCensor;thicksL];
        else
            if dvars_TM
                TempMask = dvars_TM;

                TempMask(TempMask>1) = 1;
                TempMask = conv(convVet,double(TempMask));% convolution with a Convoltion function
                TempMask = TempMask(1:dyn4,1);
                TempMask(TempMask>1) = 1;

                TempMaskIdxs = find(TempMask);
                Censu = zeros(dyn4,numel(TempMaskIdxs));
                for i = 1:numel(TempMaskIdxs)
                    Censu(TempMaskIdxs(i),i) = TempMask(TempMaskIdxs(i));
                end

                if ~isempty(Frgval)
                    Frgval_tmp = Frgval(:,1:end-1);
                    Frgval_tmp = [Frgval_tmp,Censu,ones(dyn4,1)];
                    Frgval = Frgval_tmp;
                else
                    Frgval_tmp = [Censu,ones(dyn4,1)];
                    Frgval = Frgval_tmp;
                end
                thicksLCensor(1:size(Censu,2),1) = {'Censor'};
                thicksL = [thicksLCensor;thicksL];
            else
                Censu = [];
            end
        end
    end

    if DVAR_tsV

       MeanDvarsTMP = MeanDvars;
       MeanDvarsTMP = MeanDvarsTMP./max(MeanDvarsTMP);
       MeanDvarsTMP = MeanDvarsTMP - mean(MeanDvarsTMP);

       if ~isempty(Frgval)
            Frgval_tmp = Frgval(:,1:end-1);
            Frgval_tmp = [MeanDvarsTMP,Frgval_tmp,ones(dyn4,1)];
            Frgval = Frgval_tmp;
        else
            Frgval_tmp = [MeanDvarsTMP',ones(dyn4,1)];
            Frgval = Frgval_tmp;
       end

       thicksL = [{'DVAR TS'};thicksL];

    end

%     thicksL(size(thicksL,1)+1,1) = {'Ones'};

    imagesc(Frgval);
    title('Regression Matrix','FontSize', 14);
    set(gca,'XTick',1:size(thicksL,1))
    set(gca,'XTickLabel',thicksL)
    if size(thicksL,1) <= 10
        xticklabel_rotate([],90,[],'Fontsize',11)
    end
    if size(thicksL,1) > 10 &&  size(thicksL,1) <=20
        xticklabel_rotate([],90,[],'Fontsize',10)
    end
    if size(thicksL,1) > 20
        xticklabel_rotate([],90,[],'Fontsize',8)
    end

    drawnow
    imgRR = getframe(figu);
    imwrite(imgRR.cdata, [pathFunc,dirname,filesep,filesep,'Regression Matrix.tif']);
    
end
