


% function DELTA_REVISED_n



SIG = 20;
BIN = 50;

clear SS_P_DEL_OT PVAL_DEL_OT SS_P_DEL_LL PVAL_DEL_LL

for ii=1:2
    
    
    if ii==1
        START_R =  300;   END_R = 700;
        START_T = -1400;   END_T = 200;
        START_M = -200;   END_M = 300;
    end
    
    if ii==2
        START_R =  300;   END_R = 700;
        START_T = -1400;   END_T = 400;
        START_M = -400;   END_M = 300;
    end
    
    
    n=[round(length(START_R:END_R),-2) round(length(START_T:END_T),-2) round(length(START_M:END_M),-2)]/100;
    
    TIME_R = START_R:END_R;
    TIME_T = START_T:END_T;
    TIME_M = START_M:END_M;
    
    P_R = NaN(1,length(TIME_R));
    P_T = NaN(1,length(TIME_T));
    P_M = NaN(1,length(TIME_M));
    
    
    clear allP1 allP2 allP ;
    
    if ii==1 TRIALS = 2:CHANGE;  STR = 'OT'; end
    if ii==2 TRIALS = CHANGE+1:CHANGE+29; STR = 'LL'; end
    
    SPIKES_SS = Spikes.S(TRIALS);
    
    MOVT1 = Infos(TRIALS-1,11);
    STIM2 = Infos(TRIALS,4);
    MOVT2 = Infos(TRIALS,11);
    
    PSTH_R = PSTH_RETURN_n(SPIKES_SS,MOVT1,START_R,END_R,SIG);
    PSTH_T = PSTH_RETURN_n(SPIKES_SS,STIM2,START_T,END_T,SIG);
    PSTH_M = PSTH_RETURN_n(SPIKES_SS,MOVT2,START_M,END_M,SIG);
    
    clear PSTH TIME;
    
    PSTH = [PSTH_R PSTH_T PSTH_M];
    TIME = [TIME_R TIME_T TIME_M];
    
    if ii==1 SS_PSTH_OT = nanmean(PSTH); SS_TIME_OT = (TIME); end;
    if ii==2 SS_PSTH_LL = nanmean(PSTH); SS_TIME_LL = (TIME); end;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F(ii) = figure();
    
    s2 =  subplot(4,sum(n),[1:n(1)]+sum(n)*2); s2Pos = get(s2,'position');
    S1(ii) = s2;
    hold on;
    PSTHe_n(SPIKES_SS,MOVT1,START_R,END_R,SIG,[0 0 0],1,0);
    YLIM = ylim; ylim([0 YLIM(2)]);
    hold on;
    plot([300 300],[0 0],'ok','markersize',4); xtickangle(45); 
    
    s5 =  subplot(4,sum(n),[n(1)+1:n(1)+n(2)]+sum(n)*2); s5Pos = get(s5,'position');
    S2(ii) = s5 ;
    hold on;
    PSTHe_n(SPIKES_SS,STIM2,START_T,END_T,SIG,[0 0 0],1,0);
    YLIM = ylim; ylim([0 YLIM(2)]);
    hold on;
    plot([0 0],[0 0],'ok','markersize',4)
    set(gca,'ytick',[]); set(gca,'yticklabel',[]); ylabel([]); xtickangle(45);
    
    s8 =  subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]+sum(n)*2); s8Pos = get(s8,'position');
    S3(ii) = s8 ;
    hold on;
    PSTHe_n(SPIKES_SS,MOVT2,START_M,END_M,SIG,[0 0 0],1,0);
    YLIM = ylim; ylim([0 YLIM(2)]);
    hold on;
    plot([0 0],[0 0],'ok','markersize',4);
    set(gca,'ytick',[]); set(gca,'yticklabel',[]); ylabel([]); xtickangle(45);
    plot([0 0],ylim,'-k');
    
    
    
    suptitle(STR);
    linkaxes([s8, s2, s5 ],'y');
    set(s2,'position',s2Pos); set(s5,'position',s5Pos); set(s8,'position',s8Pos);
    
    subplot(4,sum(n),[n(1)+1:n(1)+n(2)]+sum(n)*2); plot([0 0],ylim,'-k'); set(gca,'ycolor',[1 1 1])
    subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]+sum(n)*2); plot([0 0],ylim,'-k'); set(gca,'ycolor',[1 1 1])
    
    
    
    %%%%%% plot SS--------------------
    %     TRIALS = CHANGE:CHANGE+30;
    INDC = find(Infos(TRIALS,10)==1)+TRIALS(1)-1;
    INDW = find(Infos(TRIALS,10)==2)+TRIALS(1)-1;
    
    
    ss1 = subplot(4,sum(n),[1:n(1)]); ss1Pos = get(ss1,'position');
    hold on;
    PSTHe_n(Spikes.S(INDC),Infos(INDC,11),START_R,END_R,30,[0 0 1],0.6,0);
    PSTHe_n(Spikes.S(INDW),Infos(INDW,11),START_R,END_R,30,[1 0 0],0.6,0);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]);
    clear P1 P2
    P1 = PSTH_RETURN_n(Spikes.S(INDC,1),Infos(INDC,11),START_R,END_R,30);
    P2 = PSTH_RETURN_n(Spikes.S(INDW),Infos(INDW,11),START_R,END_R,30);
    allP1 = nanmean(P1); allP2 = nanmean(P2);
    parfor kk=1:size(P1,2)
        try
        P_R(1,kk) = ttest_NN(P1(:,kk),P2(:,kk));
        catch
        end;
    end
    allP = P_R;
    
    
    sx1 = subplot(4,sum(n),[1:n(1)]+sum(n)); sx1Pos = get(sx1,'position'); %%%%%%%%%%%%%% difffffffff
    hold on; 
    plot(nanmean(P1)-nanmean(P2));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]);
    xlim([1 length(P1)]);
    plot(xlim,[0 0],'-k');
    
    
    ss2 =  subplot(4,sum(n),[n(1)+1:n(1)+n(2)]); ss2Pos = get(ss2,'position');
    hold on;
    PSTHe_n(Spikes.S(INDC+1),Infos(INDC+1,4),START_T,END_T,30,[0 0 1],0.6,0);
    PSTHe_n(Spikes.S(INDW+1),Infos(INDW+1,4),START_T,END_T,30,[1 0 0],0.6,0);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); ylabel([]); set(gca,'ycolor',[1 1 1]);
    clear P1 P2
    P1 = PSTH_RETURN_n(Spikes.S(INDC+1,1),Infos(INDC+1,4),START_T,END_T,30);
    P2 = PSTH_RETURN_n(Spikes.S(INDW+1),Infos(INDW+1,4),START_T,END_T,30);
    allP1 = [allP1 nanmean(P1)]; allP2 = [allP2 nanmean(P2)];
    parfor kk=1:size(P1,2)
        try
        P_T(1,kk) = ttest_NN(P1(:,kk),P2(:,kk));
        catch
        end;
    end
    allP = [allP P_T];
    
    
    sx2 = subplot(4,sum(n),[n(1)+1:n(1)+n(2)]+sum(n)); sx2Pos = get(sx2,'position'); %%%%%%%%%%%%%% difffffffff
    hold on; 
    plot(nanmean(P1)-nanmean(P2));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]); set(gca,'yticklabel',[]); ylabel([]); set(gca,'ycolor',[1 1 1]);
    xlim([1 length(P1)]);
    plot(xlim,[0 0],'-k'); ylim([-50 50]);
    
    
    ss3 = subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]); ss3Pos = get(ss3,'position');
    hold on;
    PSTHe_n(Spikes.S(INDC+1),Infos(INDC+1,11),START_M,END_M,30,[0 0 1],0.6,0);
    PSTHe_n(Spikes.S(INDW+1),Infos(INDW+1,11),START_M,END_M,30,[1 0 0],0.6,0);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); ylabel([]);
    set(gca,'ycolor',[1 1 1]);
    clear P1 P2
    P1 = PSTH_RETURN_n(Spikes.S(INDC+1,1),Infos(INDC+1,11),START_M,END_M,30);
    P2 = PSTH_RETURN_n(Spikes.S(INDW+1),Infos(INDW+1,11),START_M,END_M,30);
    allP1 = [allP1 nanmean(P1)]; allP2 = [allP2 nanmean(P2)];
    parfor kk=1:size(P1,2)
        try
        P_M(1,kk) = ttest_NN(P1(:,kk),P2(:,kk));
        catch
        end;
    end
    allP = [allP P_M];
    
    
    
    sx3 = subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]+sum(n)); sx3Pos = get(sx3,'position'); %%%%%%%%%%%%%% difffffffff
    hold on; 
    plot(nanmean(P1)-nanmean(P2));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); xlabel([]); set(gca,'yticklabel',[]); ylabel([]); set(gca,'ycolor',[1 1 1]);
    xlim([1 length(P1)]);
    plot(xlim,[0 0],'-k');
    
    
    
    subplot(4,sum(n),[1:n(1)]); ylim([round(nanmax(nanmin([allP1 allP2])-5,0)) round(nanmax([allP1 allP2])+5)]);
    %     yticks([round(nanmax(nanmin([allP1 allP2])-5,0)) round(nanmax([allP1 allP2])+5)]);
    
    subplot(4,sum(n),[n(1)+1:n(1)+n(2)]); ylim([round(nanmax(nanmin([allP1 allP2])-5,0)) round(nanmax([allP1 allP2])+5)]);
    subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]); ylim([round(nanmax(nanmin([allP1 allP2])-5,0)) round(nanmax([allP1 allP2])+5)]);
    
    linkaxes([ ss3, ss1, ss2],'y'); linkaxes([ sx3, sx1, sx2],'y');
    set(ss1,'position',ss1Pos); set(ss2,'position',ss2Pos); set(ss3,'position',ss3Pos);
    YLIM = ylim;
    subplot(4,sum(n),[1:n(1)]); hold on; clear ind; ind = find(P_R<0.005); plot(TIME_R(ind),YLIM(2)*ones(size(ind)),'.k');
    subplot(4,sum(n),[n(1)+1:n(1)+n(2)]); hold on; plot([0 0],ylim,'-k'); clear ind; ind = find(P_T<0.005); plot(TIME_T(ind),YLIM(2)*ones(size(ind)),'.k');
    subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]); hold on; plot([0 0],ylim,'-k'); clear ind; ind = find(P_M<0.005); plot(TIME_M(ind),YLIM(2)*ones(size(ind)),'.k');
    
    
    subplot(4,sum(n),[n(1)+n(2)+1:sum(n)]+sum(n)*2); hold on; plot([END_M-100 END_M],[0 0],'-K','linewidth',2);
    
    
    if ii==1
        SS_P_DEL_OT(1,:) = allP1;
        SS_P_DEL_OT(2,:) = allP2;
        PVAL_DEL_OT = allP;
    end
    
    if ii==2
        SS_P_DEL_LL(1,:) = allP1;
        SS_P_DEL_LL(2,:) = allP2;
        PVAL_DEL_LL = allP;
    end
    
    clear allP1 allP2 allP
    
end


linkaxes([S1(1), S2(1), S3(1),S1(2), S2(2), S3(2)],'y');



cd(Results_dir)

filename = strcat('SS_epoch_OT'); print(F(1), '-dpdf', filename, '-r600');
filename = strcat('SS_epoch_LL'); print(F(2), '-dpdf', filename, '-r600');
clear F;


newDEL.SS_P_DEL_OT = SS_P_DEL_OT;
newDEL.PVAL_DEL_OT = PVAL_DEL_OT;
newDEL.SS_P_DEL_LL = SS_P_DEL_LL;
newDEL.PVAL_DEL_LL = PVAL_DEL_LL;



save(POP_file,'SS_PSTH_OT','SS_PSTH_LL','SS_TIME_OT','SS_TIME_LL','newDEL','-append');



% end
