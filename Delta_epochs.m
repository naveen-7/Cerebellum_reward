 

 %Plots and caluclates delta for each input epoch
% written by naveen at cumc on 8/12/17
 % Delta_epochs
    
    t1 = EPOCH_START; t2 = EPOCH_END;
    if strcmp(ALIGN,'T')
        P_use = P_CN_T_smooth;
        P_use_BEF = P_CN_T_smooth_BEF;
        TIME = time_T;
    elseif strcmp(ALIGN,'M')
        P_use = P_CN_M_smooth;
        P_use_BEF = P_CN_M_smooth_BEF;
        TIME = time_M;
    end
    
    
    clear VALUE_C VALUE_C_RANDOM;
    for i = 1:size(P_CN_T_smooth,1)
        VALUE_C(i,1) = nanmean(P_use(i,find(TIME==t1):find(TIME==t2))); % avg
        VALUE_C(i,2) = nanstd(P_use(i,find(TIME==t1):find(TIME==t2))); %sqrt(find(TIME==t2)-find(TIME==t1)); % stddev
        VALUE_C(i,3) = nanstd(P_use(i,find(TIME==t1):find(TIME==t2)))/sqrt(w_bin); %sqrt(find(TIME==t2)-find(TIME==t1)); % stderror
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    VALUE_C_BEF(1,1) = nanmean(P_use_BEF(1,find(TIME==t1):find(TIME==t2))); % avg
    VALUE_C_BEF(1,2) = nanstd(P_use_BEF(1,find(TIME==t1):find(TIME==t2))); % stddev
    VALUE_C_BEF(1,3) = nanstd(P_use_BEF(1,find(TIME==t1):find(TIME==t2)))/sqrt(w_bin); % std error
    

    F = figure();
    WID = 12;
    LEN = 3;
    FS=10;
    
    ax1 = subplot(WID,LEN,[1 2 4 5 7 8 10 11 13 14]);
    hold on;
    clear ylim yLim
    plot(x_trial_BEF,per_corr_BEF,'-','color',[1	0.7255	0.0588],'LineWidth',2);
    plot(x_trial_BEF,per_corr_BEF,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
    
    plot(x_trial,per_corr,'-','color',[1	0.7255	0.0588],'LineWidth',2);
    plot(x_trial,per_corr,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
    yLim = ylim;
    if yLim(1)>20
        ylim([20 yLim(2)]);
    end
    yLim = ylim;
    plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color','k');
    %     title('Learning curve');
    ylabel('% Correct');
    box off;
    ax1.FontSize=FS;
    xlim([1 size(Infos,1)]);
    
    xlim([x_trial_BEF x_trial(length(x_trial))]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    
    XVAL_TOT = [x_trial_BEF x_trial];
    
    ax2 = subplot(WID,LEN,[16 17 19 20 22 23 25 26 28 29 ]);
    hold on;
    XVAL = x_trial;
   
    
    errorline_n(x_trial_BEF,VALUE_C_BEF(:,1),VALUE_C_BEF(:,2),2,[1 0.5 0.5],0.2,1)
    errorline_n(XVAL,VALUE_C(:,1),VALUE_C(:,2),2,[1 0.5 0.5],0.2,1)
  
    
    ylabel('Diff Firing rate (Sp/s)');
    box off;
    ax2.FontSize=FS;
    plot(xlim,[0 0],'-k','linewidth',2)
    xlabel('Trail number');
    xlim([x_trial_BEF x_trial(length(x_trial))]);

    plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color','k');
    
    
    
    %%%%%%%%%%%%%STATS%%%%%%%%%%%%%%
    ax3 = subplot(WID,LEN,[3 6 9 12 15]);
    axis off;
    text(0,1,'Prev Corr')
    text(0,0.8,strcat(num2str(EPOCH_START),'to',num2str(EPOCH_END),'ms @',ALIGN));
    XVAL_CW = XVAL;
    
    
    
    VALUE_C_TOT = [VALUE_C_BEF; VALUE_C];
    per_corr_TOT = [per_corr_BEF per_corr]';
    tempyy = VALUE_C_TOT(:,1);
    clear CORR_LC_P
    [RHO_LC_PSTH p_LC_PSTH] = corr(per_corr_TOT,tempyy);
    
    
  
    
    text(0,0.5,strcat('RHO = ',num2str(RHO_LC_PSTH)));
    text(0,0.4,strcat('p = ',num2str(p_LC_PSTH)));
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VALUE_C_BEF;
    VALUE_C_DUR(1) = nanmean(VALUE_C(1:2,1));
    VALUE_C_AFT(1) = nanmean(VALUE_C(length(VALUE_C)-1:length(VALUE_C),1));
    VALUE_C_DUR(2) = nanmean(VALUE_C(1:2,3));
    VALUE_C_AFT(2) = nanmean(VALUE_C(length(VALUE_C)-1:length(VALUE_C),3));
    
    
    ax4 = subplot(WID,LEN,[18 21 24 27 30]);
    hold on;
    
    x_plot = [1 2 3];
    y_plot = [VALUE_C_BEF(1) VALUE_C_DUR(1) VALUE_C_AFT(1)];
    e_plot = [VALUE_C_BEF(2) VALUE_C_DUR(2) VALUE_C_AFT(2)];
    
    bar(x_plot,y_plot,'facecolor',[1 0.5 0.5]);
    e = errorbar(x_plot,y_plot,e_plot,'.');
    e.Color = [0 0 0];
    
    box off;
    get(ax4,'XTickLabel');
    set(ax4,'XTickLabel',{'Bef','Dur','Aft'}) %shows 1 to 11
    ylabel('Sp/s');
    
    
    suptitle(strcat(NOME(1:8),'-',NOME(10:12),'-PSTH-SS-PREV-EPOCH-DYN-',STRING{IND}));
    
    
    cd(Results_dir)
    filename = strcat(NOME,'_PSTH_SS_PREV_EPOCH_DYN_',STRING{IND});
    print(F, '-dpdf', filename, '-r600')


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if RHO_LC_PSTH>0 & p_LC_PSTH<0.05
        GROUP = 2;
    end
    
    if RHO_LC_PSTH<0 & p_LC_PSTH<0.05
        GROUP = 1;
    end
    
    if ~(RHO_LC_PSTH>0 & p_LC_PSTH<0.05) & ~(RHO_LC_PSTH<0 & p_LC_PSTH<0.05)
        GROUP=[];
    end
    
    
    if strcmp(upper(ALIGN),'T')
        if EPOCH_START<=0 & EPOCH_END<=100
            EPOCH_TYPE = 'B';
        end
        if EPOCH_START>=0 & EPOCH_END>=100
            EPOCH_TYPE = 'T';
        end
    end
    
    
    if strcmp(upper(ALIGN),'M')
        if EPOCH_START<=50 & EPOCH_END<=300
            EPOCH_TYPE = 'M';
        end
        if EPOCH_START>=100 & EPOCH_END>=200
            EPOCH_TYPE = 'P';
        end
    end
    
    GROUP
    EPOCH_TYPE
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    

 %% SAVING ZONE -------------------------

 
 
 
 DELTA.MAINcorr_RHO(IND) = RHO_LC_PSTH;
 DELTA.MAINcorr_P(IND) = p_LC_PSTH;
 
 DELTA.QUANT{IND}=[VALUE_C_BEF(1) VALUE_C_DUR(1) VALUE_C_AFT(1)];
 
 DELTA.LCx{IND} = [x_trial_BEF x_trial];
 DELTA.LCy{IND} = [per_corr_BEF per_corr];
 DELTA.delta{IND} = [VALUE_C_BEF(:,1) VALUE_C(:,1)'];

%%