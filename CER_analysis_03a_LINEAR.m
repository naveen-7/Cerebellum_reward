

%% function CER_analysis_03a_LINEAR

% Gross analysis of prev wrong and prev corr for
% before, during and after learning



%% Wrong -> Anything


Type3 = NaN(length(Infos),1);

for i=2:length(Infos)
    
    if (Infos(i-1,10)==2)
        Type3(i,1)=1;
    end
    
end




%% PRIOR WRONG TRIALS -------------------------------------------

clear PSTH_bef PSTH_aft PSTH_bef_std PSTH_aft_std
clear points_bef_T points_aft_T points_bef_R points_aft_R



for XX =1:2   % Running twice; each time align differently
    
    if XX ==1 Align_code = 4; end
    if XX ==2 Align_code = 11; end
    
    if Align_code == 11
        Start = -750; End = 750;
    end
    
    if Align_code == 4
        Start = -450; End = 1050;
    end
    
    
    
    
    % BEFORE CHANGE ----------

    clear temp_trials
    temp_trials = find(Type3(1:CHANGE,1)>=1);
    no_trials = length(temp_trials)
    no_trials_BC_W = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,60,[0.6353    0.8039    0.3529]);
%         delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    
    for j=1:size(PSTH_type1,2)
        PSTH_bef(1,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(1,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(1,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(1,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    
    NUM_TRLS(1,1) = no_trials;
    PSTH_bef_BC_SS     = PSTH_bef(1,:);
    PSTH_aft_BC_SS     = PSTH_aft(1,:);
    PSTH_bef_std_BC_SS = PSTH_bef_std(1,:);
    PSTH_aft_std_BC_SS = PSTH_aft_std(1,:);
    
    
    
    
    % DURING CHANGE ----------
   
    clear temp_trials
%     temp_trials = find(Type3(CHANGE:LEARNT,1)>=1)+CHANGE-1;
    temp_trials = find(Type3(CHANGE:CHANGE+15,1)>=1)+CHANGE-1;
    no_trials = length(temp_trials)
    no_trials_DC_W = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,60,[0.6353    0.8039    0.3529]);
%         delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    
    for j=1:size(PSTH_type1,2)
        PSTH_bef(2,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(2,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(2,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(2,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    
    NUM_TRLS(2,1) = no_trials;
    PSTH_bef_DC_SS     = PSTH_bef(2,:);
    PSTH_aft_DC_SS     = PSTH_aft(2,:);
    PSTH_bef_std_DC_SS = PSTH_bef_std(2,:);
    PSTH_aft_std_DC_SS = PSTH_aft_std(2,:);
    
    
    % AFTER CHANGE ----------
        
    clear temp_trials
    temp_trials = find(Type3(LEARNT+5:size(Infos,1),1)>=1)+LEARNT;
    no_trials = length(temp_trials)
    no_trials_AC_W = no_trials;
    
    PSTH_type1 =  NaN(no_trials,(End-Start));
    PSTH_type2 =  NaN(no_trials,(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,60,[0.9804    0.5020    0.4471]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,60,[0.6353    0.8039    0.3529]);
%         delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    for j=1:size(PSTH_type1,2)
        PSTH_bef(3,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(3,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(3,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(3,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    NUM_TRLS(3,1) = no_trials;
    PSTH_bef_AC_SS     = PSTH_bef(3,:);
    PSTH_aft_AC_SS     = PSTH_aft(3,:);
    PSTH_bef_std_AC_SS = PSTH_bef_std(3,:);
    PSTH_aft_std_AC_SS = PSTH_aft_std(3,:);
    
    
    
    
    
    % Plotting subplot
    clear fff;
    fff = figure();
    hold on;
    clear ylim yLim;
    
    
    for kk=1:3
        
        subplot(3,3,kk)
        hold on;
        time = Start:End;
        
        errorline_n(time,PSTH_bef(kk,:),PSTH_bef_std(kk,:),1,[0.9804    0.5020    0.4471],0.5)
        errorline_n(time,PSTH_aft(kk,:),PSTH_bef_std(kk,:),1,[0.2745    0.5098    0.7059],0.5)
        
        xlim([time(1)+50 time(length(time))-50]);
        xlabel('Time in ms','FontSize',12);
        ylabel('Sp/s','FontSize',12);
        
        if kk==1 title('Before Change','FontSize',12); end
        if kk==2 title('During learning','FontSize',12); end
        if kk==3 title('After learning','FontSize',12); end
        
        set(gca,'FontSize',12,'LineWidth',1)
        set(gcf, 'PaperUnits','inches','PaperSize',[8 8],'PaperPosition',[1 1 6.65 5])
        box off;
        
        yLim(kk,:) = ylim;
        
    end
    
    hold on;
    minylim = min(yLim(:,1));
    maxylim = max(yLim(:,2));
    
    for kk=1:3
        subplot(3,3,kk)
        ylim([minylim maxylim]);
%         ylim([20 60]);
        plot([0 0],ylim,'--','color',[0 0 0],'LineWidth',1);
    end
    
    
    
    
    % PLOTTING THE DIFFERENCES -----------------
    DIFF_BL = PSTH_aft_BC_SS-PSTH_bef_BC_SS;   % Difference during learning aligned to target
    DIFF_DL = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
    DIFF_AL = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    
    
    subplot(3,3,4)
    hold on;
    plot(time,DIFF_BL);
    yLim(3,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    subplot(3,3,5)
    hold on;
    plot(time,DIFF_DL);
    yLim(3,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    subplot(3,3,6)
    hold on;
    plot(time,DIFF_AL);
    yLim(4,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    
    minylim = min(yLim(3,1),yLim(4,1));
    maxylim = max(yLim(3,2),yLim(4,2));
    
    for kk=4:6
        subplot(3,3,kk)
        ylim([minylim maxylim]);
        plot([0 0],ylim,'--','color',[0 0 0],'LineWidth',1);
    end
    
    
    
    
    if Align_code == 4
        suptitle('learning SS @T');
        filename = 'WN_LEARNING_SS_@T';
        cd(Results_dir);
        print(fff, '-dpdf', filename, '-r400');
    end
    
    if Align_code == 11
        suptitle('learning SS @R');
        filename = 'WN_LEARNING_SS_@R';
        cd(Results_dir);
        print(fff, '-dpdf', filename, '-r400');
    end
    
    
    
    % POP PREP
    
    if Align_code == 4
        DIFF_BL_T_WRG = PSTH_aft_BC_SS-PSTH_bef_BC_SS;
        DIFF_DL_T_WRG = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
        DIFF_AL_T_WRG = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    end
    
    if Align_code == 11
        DIFF_BL_R_WRG = PSTH_aft_BC_SS-PSTH_bef_BC_SS;
        DIFF_DL_R_WRG = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
        DIFF_AL_R_WRG = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    end
    
    
    
    
    if Align_code == 4
        PSTH_bef_BC_SS_WRG_T = PSTH_bef_BC_SS;
        PSTH_aft_BC_SS_WRG_T = PSTH_aft_BC_SS;
        PSTH_bef_DC_SS_WRG_T = PSTH_bef_DC_SS;
        PSTH_aft_DC_SS_WRG_T = PSTH_aft_DC_SS;
        PSTH_bef_AC_SS_WRG_T = PSTH_bef_AC_SS;
        PSTH_aft_AC_SS_WRG_T = PSTH_aft_AC_SS;
    end
    
    if Align_code == 11
        PSTH_bef_BC_SS_WRG_R = PSTH_bef_BC_SS;
        PSTH_aft_BC_SS_WRG_R = PSTH_aft_BC_SS;
        PSTH_bef_DC_SS_WRG_R = PSTH_bef_DC_SS;
        PSTH_aft_DC_SS_WRG_R = PSTH_aft_DC_SS;
        PSTH_bef_AC_SS_WRG_R = PSTH_bef_AC_SS;
        PSTH_aft_AC_SS_WRG_R = PSTH_aft_AC_SS;
    end
    
    
    
end


save(MERGE_file,'DIFF_BL_T_WRG','DIFF_BL_R_WRG','DIFF_DL_T_WRG','DIFF_AL_T_WRG','DIFF_DL_R_WRG','DIFF_AL_R_WRG','-append');
save(ALLCELLS_file,'DIFF_BL_T_WRG','DIFF_BL_R_WRG','DIFF_DL_T_WRG','DIFF_AL_T_WRG','DIFF_DL_R_WRG','DIFF_AL_R_WRG','-append');






















%% PRIOR CORRECT TRIALS -------------------------------------------

Type8 = NaN(length(Infos),1);

for i=2:length(Infos)
    
    if (Infos(i-1,10)==1)
        Type8(i,1)=1;
    end
    
end



clear PSTH_bef PSTH_aft PSTH_bef_std PSTH_aft_std
clear points_bef_T points_aft_T points_bef_R points_aft_R



for XX =1:2   % Running twice; each time align differently
    
    if XX ==1 Align_code = 4; end
    if XX ==2 Align_code = 11; end
    
    if Align_code == 11
        Start = -750; End = 750;
        time_M = Start:End;
    end
    
    if Align_code == 4
        Start = -450; End = 1050;
        time_T = Start:End;
    end
    
    
    
    
    % BEFORE CHANGE ----------

    clear temp_trials
    temp_trials = find(Type8(1:CHANGE,1)>=1);
    no_trials = length(temp_trials)
    no_trials_BC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,60,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,60,[0.6353    0.8039    0.3529]);
        delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    
    for j=1:size(PSTH_type1,2)
        PSTH_bef(1,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(1,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(1,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(1,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    
    NUM_TRLS(1,1) = no_trials;
    PSTH_bef_BC_SS     = PSTH_bef(1,:);
    PSTH_aft_BC_SS     = PSTH_aft(1,:);
    PSTH_bef_std_BC_SS = PSTH_bef_std(1,:);
    PSTH_aft_std_BC_SS = PSTH_aft_std(1,:);
    
    
    
    
    % DURING CHANGE ----------
    clear temp_trials
    temp_trials = find(Type8(CHANGE:CHANGE+15,1)>=1)+CHANGE-1;
    no_trials = length(temp_trials)
    no_trials_DC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,1+(End-Start));
    PSTH_type2 =  NaN(no_trials,1+(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,40,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,40,[0.6353    0.8039    0.3529]);
        delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    

    for j=1:size(PSTH_type1,2)
        PSTH_bef(2,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(2,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(2,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(2,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    
    NUM_TRLS(2,1) = no_trials;
    PSTH_bef_DC_SS     = PSTH_bef(2,:);
    PSTH_aft_DC_SS     = PSTH_aft(2,:);
    PSTH_bef_std_DC_SS = PSTH_bef_std(2,:);
    PSTH_aft_std_DC_SS = PSTH_aft_std(2,:);
    
    
    % AFTER CHANGE ----------
        
    clear temp_trials
    temp_trials = find(Type8(LEARNT+5:size(Infos,1),1)>=1)+LEARNT;
    no_trials = length(temp_trials)
    no_trials_AC_C = no_trials;
    
    PSTH_type1 =  NaN(no_trials,(End-Start));
    PSTH_type2 =  NaN(no_trials,(End-Start));
    
    clear C_tim;
    for xy = 1:length(temp_trials)  % Running for all trials in each case
        Sig = Spikes.S{temp_trials(xy)-1,1};
        PSTH1 = PSTH_ONE_n(Sig,Infos(temp_trials(xy)-1,Align_code),Start,End,60,[0.4000    0.8039         0]);
        Sig = Spikes.S{temp_trials(xy),1};
        PSTH2 = PSTH_ONE_n(Sig,Infos(temp_trials(xy),Align_code),Start,End,60,[0.6353    0.8039    0.3529]);
        delete(gcf); delete(gcf);
        PSTH_type1(xy,1:length(PSTH1)) = PSTH1;
        PSTH_type2(xy,1:length(PSTH2)) = PSTH2;
    end
    
    for j=1:size(PSTH_type1,2)
        PSTH_bef(3,j) = nanmean(PSTH_type1(:,j));
        PSTH_aft(3,j) = nanmean(PSTH_type2(:,j));
        PSTH_bef_std(3,j) = nanstd(PSTH_type1(:,j))/sqrt(no_trials);
        PSTH_aft_std(3,j) = nanstd(PSTH_type2(:,j))/sqrt(no_trials);
    end
    
    NUM_TRLS(3,1) = no_trials;
    PSTH_bef_AC_SS     = PSTH_bef(3,:);
    PSTH_aft_AC_SS     = PSTH_aft(3,:);
    PSTH_bef_std_AC_SS = PSTH_bef_std(3,:);
    PSTH_aft_std_AC_SS = PSTH_aft_std(3,:);
    
    
    
    
    
    % Plotting subplot
    clear fff;
    fff = figure();
    hold on;
    clear ylim yLim;
    
    
    for kk=1:3
        
        subplot(3,3,kk)
        hold on;
        time = Start:End;
        errorline_n(time,PSTH_bef(kk,:),PSTH_bef_std(kk,:),1,[0.4000    0.8039         0],0.5)
        errorline_n(time,PSTH_aft(kk,:),PSTH_bef_std(kk,:),1,[0.2745    0.5098    0.7059],0.5)
       
        xlim([time(1)+50 time(length(time))-50]);
        xlabel('Time in ms','FontSize',12);
        ylabel('Sp/s','FontSize',12);
        
        if kk==1 title('Before Change','FontSize',12); end
        if kk==2 title('During learning','FontSize',12); end
        if kk==3 title('After learning','FontSize',12); end
        
        set(gca,'FontSize',12,'LineWidth',1)
        set(gcf, 'PaperUnits','inches','PaperSize',[8 8],'PaperPosition',[1 1 6.65 5])
        box off;
        
        yLim(kk,:) = ylim;
       
    end
    
    hold on;
    minylim = min(yLim(:,1));
    maxylim = max(yLim(:,2));
    
    for kk=1:3
        subplot(3,3,kk)
%         ylim([minylim maxylim]);
        ylim([40 90]);
        plot([0 0],ylim,'--','color',[0 0 0],'LineWidth',1);
    end
    
    
    
    
    % PLOTTING THE DIFFERENCES -----------------
    DIFF_BL = PSTH_aft_BC_SS-PSTH_bef_BC_SS;   % Difference during learning aligned to target
    DIFF_DL = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
    DIFF_AL = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    
    
    subplot(3,3,4)
    hold on;
    plot(time,DIFF_BL);
    yLim(3,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    subplot(3,3,5)
    hold on;
    plot(time,DIFF_DL);
    yLim(3,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    subplot(3,3,6)
    hold on;
    plot(time,DIFF_AL);
    yLim(4,:) = ylim;
    xlim([time(1)+50 time(length(time))-50]);
    
    minylim = min(yLim(3,1),yLim(4,1));
    maxylim = max(yLim(3,2),yLim(4,2));
    
    for kk=4:6
        subplot(3,3,kk)
        ylim([minylim maxylim]);
        plot([0 0],ylim,'--','color',[0 0 0],'LineWidth',1);
    end
    
   
    
    if Align_code == 4
        suptitle('learning SS @T');
        filename = 'CN_LEARNING_SS_@T';
        cd(Results_dir);
        print(fff, '-dpdf', filename, '-r400');
    end
    
    if Align_code == 11
        suptitle('learning SS @R');
        filename = 'CN_LEARNING_SS_@R';
        cd(Results_dir);
        print(fff, '-dpdf', filename, '-r400');
    end
    
    
    
    % POP PREP
    
    if Align_code == 4
        DIFF_BL_T_CORR = PSTH_aft_BC_SS-PSTH_bef_BC_SS;
        DIFF_DL_T_CORR = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
        DIFF_AL_T_CORR = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    end
    
    if Align_code == 11
        DIFF_BL_R_CORR = PSTH_aft_BC_SS-PSTH_bef_BC_SS;
        DIFF_DL_R_CORR = PSTH_aft_DC_SS-PSTH_bef_DC_SS;   % Difference during learning aligned to target
        DIFF_AL_R_CORR = PSTH_aft_AC_SS-PSTH_bef_AC_SS;
    end
    
    
    
    
    if Align_code == 4
        PSTH_bef_BC_SS_CORR_T = PSTH_bef_BC_SS;
        PSTH_aft_BC_SS_CORR_T = PSTH_aft_BC_SS;
        PSTH_bef_DC_SS_CORR_T = PSTH_bef_DC_SS;
        PSTH_aft_DC_SS_CORR_T = PSTH_aft_DC_SS;
        PSTH_bef_AC_SS_CORR_T = PSTH_bef_AC_SS;
        PSTH_aft_AC_SS_CORR_T = PSTH_aft_AC_SS;
    end
    
    if Align_code == 11
        PSTH_bef_BC_SS_CORR_R = PSTH_bef_BC_SS;
        PSTH_aft_BC_SS_CORR_R = PSTH_aft_BC_SS;
        PSTH_bef_DC_SS_CORR_R = PSTH_bef_DC_SS;
        PSTH_aft_DC_SS_CORR_R = PSTH_aft_DC_SS;
        PSTH_bef_AC_SS_CORR_R = PSTH_bef_AC_SS;
        PSTH_aft_AC_SS_CORR_R = PSTH_aft_AC_SS;
    end
    
    
    
end


save(MERGE_file,'DIFF_BL_T_CORR','DIFF_BL_R_CORR','DIFF_DL_T_CORR','DIFF_AL_T_CORR','DIFF_DL_R_CORR','DIFF_AL_R_CORR','-append');
save(ALLCELLS_file,'DIFF_BL_T_CORR','DIFF_BL_R_CORR','DIFF_DL_T_CORR','DIFF_AL_T_CORR','DIFF_DL_R_CORR','DIFF_AL_R_CORR','-append');




% end
