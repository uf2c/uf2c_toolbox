function CircConnTscoreFDR(matrixR2,limiteG,ScreSize,cc,g1_VAR2,g2_VAR2,...
                    nEle,LoopAte,fact,str_nEle,Pva,tarNetLog,TarNetC,...
                    DirF,NetLabels,crit_p,g1_mean,g2_mean)
          
    RAIO = limiteG/2.2;
    ang=0:0.01:2*pi; 
    xp=RAIO*cos(ang);
    yp=RAIO*sin(ang);
    xpt=(limiteG/2.15)*cos(ang);
    ypt=(limiteG/2.15)*sin(ang);

    AnoRa = size(xp,2)/size(matrixR2,2);
    xp2 = xp(1:AnoRa:end);
    yp2 = yp(1:AnoRa:end);

    xp2t = xpt(1:AnoRa:end);
    yp2t = ypt(1:AnoRa:end);
    xp2t = xp2t+(limiteG/2);
    yp2t = yp2t+(limiteG/2);
    xp2t = xp2t(end:-1:1);
    yp2t = yp2t(end:-1:1);

    xp2 = xp2+(limiteG/2);
    yp2 = yp2+(limiteG/2);
    yp2 = yp2(end:-1:1);
    xp2 = xp2(end:-1:1);

    minH = (pdist([xp2(1),yp2(1);xp2(2),yp2(2)],'euclidean'));           
    Diam = max(xp) - min(xp);

    circTS2 = figure('Visible','on');
        set(circTS2,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                  'Color',[0 0 0],'Units','Normalized');
    set(gca,'xticklabels',[])
    set(gca,'yticklabels',[])
    axis off
    hold on

    xlim([-1 limiteG+1]);
    ylim([-1 limiteG+1]);

    pos1 = [0.15 0.15 0.7 0.7];
    subplot('Position',pos1,'Color',[0 0 0]);
    daspect([1 1 1])
    cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);
    axis off

     mTextBox1 = uicontrol('style','text',...
         'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
         'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
         'ForegroundColor',[1 .949 .867],'FontWeight','bold');
     set(mTextBox1,'String','Circular Connectome - DRAWING........');

    TvaCrit = tinv(1-(crit_p*.5),(size(g1_VAR2,3)+size(g2_VAR2,3)-2));
    crit_pTxt = (round(crit_p*10000))/10000;

     mTextBox2 = uicontrol('style','text',...
         'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.4  ScreSize(1)*.0125]),...
         'FontSize',15,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
         'ForegroundColor',[1 .949 .867],'FontWeight','bold');
     set(mTextBox2,'String',['Corrected results (FDR alpha ' Pva,...
         ', Uncorr. critical p: ',num2str(crit_pTxt),', |T-score|>',num2str(TvaCrit),')']);
    drawnow

    clear Tva TvaMax

    ccJet = jet(64);
    ccGre = gray(64);
    for xcv = 1:size(g1_VAR2,1)
        for ycv = 1:size(g1_VAR2,2)
            Tva(xcv,ycv) = tinv(1-(matrixR2(xcv,ycv)*.5),(size(g1_VAR2,3)+size(g2_VAR2,3)-2));
        end
    end
    Tva(Tva==inf) = 0;
    TvaMax = max(max(Tva));

    cbm1 = 0;
    cbm2 = 0;
    cbm3 = 0;

    for i = fact:LoopAte
        drawnow
        tester = i<=(str_nEle+nEle)-1;
        networkNumb = find(tester==1,1);
        colorV = cc(networkNumb,:);
        clear vet
        if i<10
            str = ['0' num2str(i)];
        else
            str = num2str(i);
        end

        vet = matrixR2(i-(fact-1),:);
        vet = find(vet);

        if ~isempty(vet)
            bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
            Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
            if abs(Angle)>pi/2
                Angle = Angle - pi;
                set(bt,'HorizontalAlignment','right')
            end
            set(bt,'Rotation',rad2deg(Angle))
            for j = 1:size(vet,2)

                if (g1_mean(i-(fact-1),vet(j))*g2_mean(i-(fact-1),vet(j)))>0
                    if abs(g1_mean(i-(fact-1),vet(j)))>abs(g2_mean(i-(fact-1),vet(j)))
                        colorVlidx = (round((Tva(i-(fact-1),vet(j))/TvaMax)*32))+32;
                        colorVl = ccJet(colorVlidx,:);
                    else
                        colorVlidx = round(((-Tva(i-(fact-1),vet(j))/TvaMax)*32))+32;
                        colorVl = ccJet(colorVlidx,:);
                    end
                else
                    colorVlidx = round((Tva(i-(fact-1),vet(j))/TvaMax)*64);
                    colorVl = ccGre(colorVlidx,:);
                end

                P1 = [xp2(i),yp2(i)];
                P2 = [xp2(vet(j)),yp2(vet(j))];

                ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                R = ((ARCalpha^2) + (4*(h^2)))/(8*h);

                hold on
                NewC = circ_cent(P1,P2,R);
                Marg = (limiteG-2*RAIO)/2;
                DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                if DisC1>DisC2
                    centerIdx = 1;
                else
                    centerIdx = 2;
                end

                AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                if sum(abs(AngRad))>pi && any(AngRad<0)
                    idxCic = find(AngRad<0);
                    AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                end

                [AngRad,ksord] = sort(AngRad);

                if isequal(ksord,[1 2])
                    TranOrd = 'descend';
                else
                    TranOrd = 'ascend';
                end

                ang2 = AngRad(1):0.03:AngRad(2);
                ang2 = [ang2,AngRad(2)];

                xunit = R * cos(ang2) + NewC(centerIdx,1);
                yunit = R * sin(ang2) + NewC(centerIdx,2);
                
                figure(circTS2)
                h1 = plot(xunit,yunit,'Color',colorVl,'LineWidth',3);
                h1.Color(4) = 0.4;

                if tarNetLog
                    testerj = vet(j)<=(str_nEle+nEle)-1;
                    networkNumbj = find(testerj==1,1);
                    colorVj = cc(networkNumbj,:);
                    if networkNumbj ~= TarNetC
                        if vet(j)<10
                            strf = ['0' num2str(vet(j))];
                        else
                            strf = num2str(vet(j));
                        end
                        bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                            colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                            'VerticalAlignment','middle');

                        Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                        if abs(Angle)>pi/2
                            Angle = Angle - pi;
                            set(bt,'HorizontalAlignment','right')
                        end
                        set(bt,'Rotation',rad2deg(Angle))
                    end
                end
            end
        end
    end

    clear A_cell Labx1
    subplot('Position',[0.04 0.01 0.9 0.015])
    imagesc(-TvaMax:(TvaMax/32):TvaMax);
    colormap(gca,'jet');
    set(gca, 'XAxisLocation', 'top','YColor',[0 0 0],'XColor',[1 1 1])
    xticks([1:64/16:64]);
    Labx1 = [(round(10*(-TvaMax:(TvaMax/8):TvaMax)))/10 TvaMax];
    A_cell = split(cellstr(num2str(Labx1)));
    xticklabels(A_cell)
    yticklabels('')
    title(gca,'T-Score G1>G2 (hot colors) G2>G1 (cool colors)','Color','w')

    set(mTextBox1,'String','Circular Connectome');
    drawnow

    set(circTS2, 'InvertHardCopy', 'off');
    imgRR = getframe(circTS2);
    imwrite(imgRR.cdata, [DirF filesep 'FDR_corrected_T-Score_Circ_Map' '.png']);
    clear imgRR
    close(circTS2)

end
