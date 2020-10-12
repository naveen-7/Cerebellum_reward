%% CER_analysis_03b_LINEAR

% DIfference curing with learning

% load(ALLCELLS_file);


Sigma=50;
Signal = Spikes.S;

for I=1:2
    if I==1 Align_code = 4; end
    if I==2 Align_code = 11; end
    
    if Align_code == 11
        Start = -750; End = 750;
        time_M = Start:End;
    end
    
    if Align_code == 4
        Start = -450; End = 1050;
        time_T = Start:End;
    end
    
    clear P
    for i=1:size(Infos,1)
        P(i,:) = PSTH_ONE_n(Signal{i,1},Infos(i,Align_code),Start,End,Sigma);
    end
    
    if I==1 P_T = P; end
    if I==2 P_M = P; end
end


% Take bins from CHANGE to end of recording
% Take the diff between curr and prior wrong ||ly for prior corr and plot


n_bin = 0;  % number of bins
w_bin = round(0.1*(size(Infos,1))); % bin width
s_bin = round(0.05*(size(Infos,1)));  % shift in bin

% % w_bin = 10; % bin width
% % s_bin = 5;  % shift in bin

clear IND_wn IND_cn P_WN_CURR_T P_WN_PREV_T P_WN_CURR_M P_WN_PREV_M ...
    P_CN_CURR_T P_CN_PREV_T P_CN_CURR_M P_CN_PREV_M P_CN_T P_WN_T P_CN_M P_WN_M

clear per_corr x_trial RT

for i=CHANGE:s_bin:size(Infos,1)
    if i <= size(Infos,1)+1-w_bin
        i
        
        n_bin = n_bin+1;
        
        count = 0;
        for j=i:i+w_bin-1
            if Infos(j,10)==1
                count = count+1;
            end
        end
        per_corr(n_bin) = (count/w_bin)*100;
        x_trial(n_bin) = i;
        RT(n_bin) = nanmean(Infos(i:i+w_bin-1,14));
        
        
        IND_wn = find( (Infos(i:i+w_bin-1,16)==1) | (Infos(i:i+w_bin-1,16)==4))+i-1; % curr trials whose prior tials were wrong
        IND_cn = find( (Infos(i:i+w_bin-1,16)==2) | (Infos(i:i+w_bin-1,16)==3))+i-1; % curr trials whose prior tials were corr
        
        %         if length(IND_wn)>1
        %             P_WN_CURR_T(n_bin,:) = nanmean(P_T(IND_wn,:));
        %             P_WN_PREV_T(n_bin,:) = nanmean(P_T(IND_wn-1,:));
        %             P_WN_CURR_M(n_bin,:) = nanmean(P_M(IND_wn,:));
        %             P_WN_PREV_M(n_bin,:) = nanmean(P_M(IND_wn-1,:));
        %         elseif length(IND_wn)==1
        %             P_WN_CURR_T(n_bin,:) = (P_T(IND_wn,:));
        %             P_WN_PREV_T(n_bin,:) = (P_T(IND_wn-1,:));
        %             P_WN_CURR_M(n_bin,:) = (P_M(IND_wn,:));
        %             P_WN_PREV_M(n_bin,:) = (P_M(IND_wn-1,:));
        %         elseif length(IND_wn)==0
        %             P_WN_CURR_T(n_bin,:) = NaN;
        %             P_WN_PREV_T(n_bin,:) = NaN;
        %             P_WN_CURR_M(n_bin,:) = NaN;
        %             P_WN_PREV_M(n_bin,:) = NaN;
        %         end
        
        
        if length(IND_cn)>1
            P_CN_CURR_T(n_bin,:) = nanmean(P_T(IND_cn,:));
            P_CN_PREV_T(n_bin,:) = nanmean(P_T(IND_cn-1,:));
            P_CN_CURR_M(n_bin,:) = nanmean(P_M(IND_cn,:));
            P_CN_PREV_M(n_bin,:) = nanmean(P_M(IND_cn-1,:));
        elseif length(IND_cn)==1
            P_CN_CURR_T(n_bin,:) = (P_T(IND_cn,:));
            P_CN_PREV_T(n_bin,:) = (P_T(IND_cn-1,:));
            P_CN_CURR_M(n_bin,:) = (P_M(IND_cn,:));
            P_CN_PREV_M(n_bin,:) = (P_M(IND_cn-1,:));
        elseif length(IND_cn)==0
            P_CN_CURR_T(n_bin,:) = NaN;
            P_CN_PREV_T(n_bin,:) = NaN;
            P_CN_CURR_M(n_bin,:) = NaN;
            P_CN_PREV_M(n_bin,:) = NaN;
        end
        
        
        % calculating the differences ------------------------------
        
        % %         P_WN_T(n_bin,:) = P_WN_CURR_T(n_bin,:)-P_WN_PREV_T(n_bin,:);
        % %         P_WN_M(n_bin,:) = P_WN_CURR_M(n_bin,:)-P_WN_PREV_M(n_bin,:);
        
        P_CN_T(n_bin,:) = P_CN_CURR_T(n_bin,:)-P_CN_PREV_T(n_bin,:);
        P_CN_M(n_bin,:) = P_CN_CURR_M(n_bin,:)-P_CN_PREV_M(n_bin,:);
        
        
    end
end



for i=1:size(P_CN_T,1)
    P_CN_T_smooth(i,:) = smooth(P_CN_T(i,:),0.3,'rloess');
    P_CN_M_smooth(i,:) = smooth(P_CN_M(i,:),0.3,'rloess');
end

% for i=1:size(P_WN_T,1)
%     if ~isnan(P_WN_T(i,:))
%         P_WN_T_smooth(i,:) = smooth(P_WN_T(i,:),0.3,'rloess');
%         P_WN_M_smooth(i,:) = smooth(P_WN_M(i,:),0.3,'rloess');
%     end
% end





%%%%%%%%%%%%%%%%%%%%
% BEFORE CHANGE ----
%%%%%%%%%%%%%%%%%%%%

count = 0;
x_trial_BEF = (CHANGE-11+CHANGE-1)/2;
RT_BEF = nanmean(Infos(CHANGE-11:CHANGE-1,14));

for i=CHANGE-11:CHANGE-1
    if Infos(i,10)==1
        count = count+1;
    end
    per_corr_BEF = (count/11)*100;
end


IND_wn = find( (Infos(CHANGE-11:CHANGE-1,16)==1) | (Infos(CHANGE-11:CHANGE-1,16)==4))+CHANGE-11; % curr trials whose prior tials were wrong
IND_cn = find( (Infos(CHANGE-11:CHANGE-1,16)==2) | (Infos(CHANGE-11:CHANGE-1,16)==3))+CHANGE-11; % curr trials whose prior tials were corr

% if length(IND_wn)>1
%     P_WN_CURR_T_BEF(1,:) = nanmean(P_T(IND_wn,:));
%     P_WN_PREV_T_BEF(1,:) = nanmean(P_T(IND_wn-1,:));
%     P_WN_CURR_M_BEF(1,:) = nanmean(P_M(IND_wn,:));
%     P_WN_PREV_M_BEF(1,:) = nanmean(P_M(IND_wn-1,:));
% elseif length(IND_wn)==1
%     P_WN_CURR_T_BEF(1,:) = (P_T(IND_wn,:));
%     P_WN_PREV_T_BEF(1,:) = (P_T(IND_wn-1,:));
%     P_WN_CURR_M_BEF(1,:) = (P_M(IND_wn,:));
%     P_WN_PREV_M_BEF(1,:) = (P_M(IND_wn-1,:));
% elseif length(IND_wn)==0
%     P_WN_CURR_T_BEF(1,:) = NaN;
%     P_WN_PREV_T_BEF(1,:) = NaN;
%     P_WN_CURR_M_BEF(1,:) = NaN;
%     P_WN_PREV_M_BEF(1,:) = NaN;
% end


if length(IND_cn)>1
    P_CN_CURR_T_BEF(1,:) = nanmean(P_T(IND_cn,:));
    P_CN_PREV_T_BEF(1,:) = nanmean(P_T(IND_cn-1,:));
    P_CN_CURR_M_BEF(1,:) = nanmean(P_M(IND_cn,:));
    P_CN_PREV_M_BEF(1,:) = nanmean(P_M(IND_cn-1,:));
elseif length(IND_cn)==1
    P_CN_CURR_T_BEF(1,:) = (P_T(IND_cn,:));
    P_CN_PREV_T_BEF(1,:) = (P_T(IND_cn-1,:));
    P_CN_CURR_M_BEF(1,:) = (P_M(IND_cn,:));
    P_CN_PREV_M_BEF(1,:) = (P_M(IND_cn-1,:));
elseif length(IND_cn)==0
    P_CN_CURR_T_BEF(1,:) = NaN;
    P_CN_PREV_T_BEF(1,:) = NaN;
    P_CN_CURR_M_BEF(1,:) = NaN;
    P_CN_PREV_M_BEF(1,:) = NaN;
end


% calculating the differences ------------------------------

% P_WN_T_BEF(1,:) = P_WN_CURR_T_BEF(1,:)-P_WN_PREV_T_BEF(1,:);
% P_WN_M_BEF(1,:) = P_WN_CURR_M_BEF(1,:)-P_WN_PREV_M_BEF(1,:);

P_CN_T_BEF(1,:) = P_CN_CURR_T_BEF(1,:)-P_CN_PREV_T_BEF(1,:);
P_CN_M_BEF(1,:) = P_CN_CURR_M_BEF(1,:)-P_CN_PREV_M_BEF(1,:);



P_CN_T_smooth_BEF(1,:) = smooth(P_CN_T_BEF(1,:),0.3,'rloess');
P_CN_M_smooth_BEF(1,:) = smooth(P_CN_M_BEF(1,:),0.3,'rloess');
%
% if ~isnan(P_WN_T_BEF(1,:))
%     P_WN_T_smooth_BEF(1,:) = smooth(P_WN_T_BEF(1,:),0.3,'rloess');
%     P_WN_M_smooth_BEF(1,:) = smooth(P_WN_M_BEF(1,:),0.3,'rloess');
% end




%%%%%%%%%%%%%%%%%%









% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % xx = P_CN_T(i,:);
% % % % % % % % % % % % % % % % xx = smooth(P_CN_T(i,:),0.3,'rloess');
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % plot(xlim,[0 0],'-k','linewidth',2)
% % % % % % % % % % % % % % % % xlim([-500 1000])
% % % % % % % % % % % % % % % title('@T');
% % % % % % % % % % % % % % % for i=1:size(P_CN_T_smooth,1)
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % %     plot(time_T,P_CN_T_smooth(i,:))
% % % % % % % % % % % % % % %     plot(xlim,[0 0],'-k','linewidth',2)
% % % % % % % % % % % % % % %     pause(1);
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % plot(xlim,[0 0],'-k','linewidth',2)
% % % % % % % % % % % % % % % % xlim([-500 1000])
% % % % % % % % % % % % % % % title('@M');
% % % % % % % % % % % % % % % for i=1:size(P_CN_M_smooth,1)
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % %     plot(time_M,P_CN_M_smooth(i,:))
% % % % % % % % % % % % % % %     plot(xlim,[0 0],'-k','linewidth',2)
% % % % % % % % % % % % % % %     pause(1);
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % end




% doing calc for plotting--------------------------

% % % P_CN_T_smooth
% % % P_CN_M_smooth
% % %
% % % P_WN_T_smooth
% % % P_WN_M_smooth


for i = 1:size(P_CN_T_smooth,1)
    
    t1 = -300; t2 = 0;
    Base_C_val(i,1) = nanmean(P_CN_T_smooth(i,find(time_T==t1):find(time_T==t2))); % avg
    Base_C_val(i,2) = nanstd(P_CN_T_smooth(i,find(time_T==t1):find(time_T==t2)))/sqrt(find(time_T==t2)-find(time_T==t1)); % std
    
    t1 = 100; t2 = 200;
    Targ_C_val(i,1) = nanmean(P_CN_T_smooth(i,find(time_T==t1):find(time_T==t2))); % avg
    Targ_C_val(i,2) = nanstd(P_CN_T_smooth(i,find(time_T==t1):find(time_T==t2)))/sqrt(find(time_T==t2)-find(time_T==t1)); % std
    
    t1 = -200; t2 = 100;
    Movt_C_val(i,1) = nanmean(P_CN_M_smooth(i,find(time_M==t1):find(time_M==t2))); % avg
    Movt_C_val(i,2) = nanstd(P_CN_M_smooth(i,find(time_M==t1):find(time_M==t2)))/sqrt(find(time_M==t2)-find(time_M==t1)); % std
    
    t1 = 200; t2 = 400;
    Pmov_C_val(i,1) = nanmean(P_CN_M_smooth(i,find(time_M==t1):find(time_M==t2))); % avg
    Pmov_C_val(i,2) = nanstd(P_CN_M_smooth(i,find(time_M==t1):find(time_M==t2)))/sqrt(find(time_M==t2)-find(time_M==t1)); % std
    
end
%
%
% if exist('P_WN_T_smooth')
%     for i = 1:size(P_WN_T_smooth,1)
%
%         t1 = -300; t2 = 0;
%         Base_W_val(i,1) = nanmean(P_WN_T_smooth(i,find(time_T==t1):find(time_T==t2))); % avg
%         Base_W_val(i,2) = nanstd(P_WN_T_smooth(i,find(time_T==t1):find(time_T==t2)))/sqrt(find(time_T==t2)-find(time_T==t1)); % std
%
%         t1 = 100; t2 = 200;
%         Targ_W_val(i,1) = nanmean(P_WN_T_smooth(i,find(time_T==t1):find(time_T==t2))); % avg
%         Targ_W_val(i,2) = nanstd(P_WN_T_smooth(i,find(time_T==t1):find(time_T==t2)))/sqrt(find(time_T==t2)-find(time_T==t1)); % std
%
%         t1 = -200; t2 = 100;
%         Movt_W_val(i,1) = nanmean(P_WN_M_smooth(i,find(time_M==t1):find(time_M==t2))); % avg
%         Movt_W_val(i,2) = nanstd(P_WN_M_smooth(i,find(time_M==t1):find(time_M==t2)))/sqrt(find(time_M==t2)-find(time_M==t1)); % std
%
%         t1 = 200; t2 = 400;
%         Pmov_W_val(i,1) = nanmean(P_WN_M_smooth(i,find(time_M==t1):find(time_M==t2))); % avg
%         Pmov_W_val(i,2) = nanstd(P_WN_M_smooth(i,find(time_M==t1):find(time_M==t2)))/sqrt(find(time_M==t2)-find(time_M==t1)); % std
%
%     end
% end




% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% % F = figure();
% % WID = 12;
% % LEN = 3;
% % FS=10;
% %
% % ax1 = subplot(WID,LEN,[1 2 4 5 7 8]);
% % hold on;
% % clear ylim yLim
% % plot(x_trial,per_corr,'-','color',[1	0.7255	0.0588],'LineWidth',2);
% % plot(x_trial,per_corr,'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% % xlabel('Trial number');
% % yLim = ylim;
% % if yLim(1)>20
% %     ylim([20 yLim(2)]);
% % end
% % yLim = ylim;
% % % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color','k');
% % title('Learning curve');
% % ylabel('% Correct');
% % box off;
% % ax1.FontSize=FS;
% % xlim([1 size(Infos,1)]);
% %
% % xlim([x_trial(1) x_trial(length(x_trial))]);
% %
% %
% %
% %
% %
% % % prev corr
% % ax2 = subplot(WID,LEN,[13 14 16 17 19 20 22 23]);
% % hold on;
% % XVAL = x_trial;
% % BASE = Base_C_val;
% % TARG = Targ_C_val;
% % MOVT = Movt_C_val;
% % PMOV = Pmov_C_val;
% %
% % MAX = nanmax(nanmax([BASE TARG MOVT PMOV]));
% % MIN = nanmin(nanmin([BASE TARG MOVT PMOV]));
% %
% % plot(xlim,[0 0],'-k','linewidth',2)
% %
% % plot(XVAL,BASE(:,1),'-','color',[0.5 0.5 0.5],'LineWidth',2);
% % plot(XVAL,BASE(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% %
% % plot(XVAL,TARG(:,1),'-','color',[0.1098    0.5255    0.9333],'LineWidth',2);
% % plot(XVAL,TARG(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.1098    0.5255    0.9333],'color','k');
% %
% % plot(XVAL,MOVT(:,1),'-','color',[0.6039    0.8039    0.1961],'LineWidth',2);
% % plot(XVAL,MOVT(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.6039    0.8039    0.1961],'color','k');
% %
% % plot(XVAL,PMOV(:,1),'-','color',[0.9333    0.4627         0],'LineWidth',2);
% % plot(XVAL,PMOV(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.9333    0.4627         0],'color','k');
% %
% % xlim([1 size(Infos,1)]);
% % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color','k');
% % ylabel('Firing rate (Sp/s)');
% % box off;
% % ax2.FontSize=FS;
% % ylim([MIN MAX])
% % plot(xlim,[0 0],'-k','linewidth',2)
% % xlabel('Trail number');
% % xlim([x_trial(1) x_trial(length(x_trial))]);
% %
% % subplot(WID,LEN,[15 18 21]);
% % axis off;
% % text(0,1,'Prev Corr')
% %
% %


% % % prev wrong
% % ax3 = subplot(WID,LEN,[25 26 28 29 31 32 34 35]);
% % hold on;
% % XVAL = x_trial(1:length(Base_W_val));
% % BASE = Base_W_val;
% % TARG = Targ_W_val;
% % MOVT = Movt_W_val;
% % PMOV = Pmov_W_val;
% %
% % MAX = nanmax(nanmax([BASE TARG MOVT PMOV]));
% % MIN = nanmin(nanmin([BASE TARG MOVT PMOV]));
% %
% % plot(xlim,[0 0],'-k','linewidth',2)
% %
% % plot(XVAL,BASE(:,1),'-','color',[0.5 0.5 0.5],'LineWidth',2);
% % plot(XVAL,BASE(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'color','k');
% %
% % plot(XVAL,TARG(:,1),'-','color',[0.1098    0.5255    0.9333],'LineWidth',2);
% % plot(XVAL,TARG(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.1098    0.5255    0.9333],'color','k');
% %
% % plot(XVAL,MOVT(:,1),'-','color',[0.6039    0.8039    0.1961],'LineWidth',2);
% % plot(XVAL,MOVT(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.6039    0.8039    0.1961],'color','k');
% %
% % plot(XVAL,PMOV(:,1),'-','color',[0.9333    0.4627         0],'LineWidth',2);
% % plot(XVAL,PMOV(:,1),'o','MarkerSize',4,'MarkerFaceColor',[0.9333    0.4627         0],'color','k');
% %
% % xlim([1 size(Infos,1)]);
% % plot([CHANGE CHANGE], ylim,'--','lineWidth',1,'color','k');
% % ylabel('Firing rate (Sp/s)');
% % box off;
% % ax4.FontSize=FS;
% % ylim([MIN MAX])
% %
% % xlabel('Trail number');
% % xlim([x_trial(1) x_trial(length(x_trial))]);
% %

% subplot(WID,LEN,[27 30 33]);
% axis off;
% text(0,1,'Prev Wrong')


% %
% % subplot(WID,LEN,[33 36]);
% % axis off;
% % text(0,0,'BASE','color',[0.5 0.5 0.5])
% % text(0,0.2,'TARG','color',[0.1098    0.5255    0.9333])
% % text(0,0.4,'MOVT','color',[0.6039    0.8039    0.1961])
% % text(0,0.6,'PMOV','color',[0.9333    0.4627         0])
% %
% %
% %
% %
% % cd(Results_dir)
% % filename = strcat(NOME,'_PSTH_SS_PREV_DYN');
% %
% % print(F, '-dpdf', filename, '-r600')







%% Epoch specific analysis --------------------

for XYXY=1:2
    
    if XYXY==1
        AAA=-450;
        BBB=1050;
        ALIGN = 'T';
        P_use = P_CN_T_smooth;
        P_use_BEF = P_CN_T_smooth_BEF;
        TIME = time_T;
    end
    
    if XYXY==2
        AAA=-750;
        BBB=750;
        ALIGN = 'M';
        P_use = P_CN_M_smooth;
        P_use_BEF = P_CN_M_smooth_BEF;
        TIME = time_M;
    end
    
    COUNT=0;
    for tt = AAA:50:BBB
        COUNT=COUNT+1;
        EPOCH_START=tt;
        EPOCH_END = EPOCH_START+200;
        t1 = EPOCH_START; t2 = EPOCH_END;
        clear VALUE_C;
        for i = 1:size(P_CN_T_smooth,1)
            VALUE_C(i,1) = nanmean(P_use(i,find(TIME==t1):find(TIME==t2))); % avg
        end
        TTEMP = nanmean(P_use_BEF(1,find(TIME==t1):find(TIME==t2))); % avg
        VALUE_C_BEF = TTEMP(1,1);
        VALUE_C_TOT = [VALUE_C_BEF; VALUE_C];
        per_corr_TOT = [per_corr_BEF per_corr]';
        tempyy = VALUE_C_TOT(:,1);
        [RHO_LC(XYXY,COUNT) p_LC(XYXY,COUNT)] = corr(per_corr_TOT,tempyy);
        
        X_TIME(XYXY,COUNT) = (EPOCH_START+EPOCH_END)/2;
    end
    
end



FFFFFF = figure();

subplot(2,2,1)
plot(X_TIME(1,:),RHO_LC(1,:))
box off;
hold on;
ylim([-1 1])
xlim([-300 1100])
ind = find(p_LC(1,:)<0.0001);
plot(X_TIME(1,ind),RHO_LC(1,ind),'*');
if isempty(ind)
    ind = find(p_LC(1,:)<0.005);
    plot(X_TIME(1,ind),RHO_LC(1,ind),'o');
end
xlabel('time from target');

subplot(2,2,2)
plot(X_TIME(2,:),RHO_LC(2,:))
box off;
hold on;
ylim([-1 1])
xlim([-600 800])
ind = find(p_LC(2,:)<0.0001);
plot(X_TIME(2,ind),RHO_LC(2,ind),'*');

if isempty(ind)
    ind = find(p_LC(2,:)<0.005);
    plot(X_TIME(2,ind),RHO_LC(2,ind),'o');
end

xlabel('time from movement');


subplot(2,2,3)
plot(X_TIME(1,:),abs(RHO_LC(1,:)))
box off;
hold on;
ylim([0 1])
xlim([-300 1100])
ind = find(p_LC(1,:)<0.0001);
plot(X_TIME(1,ind),abs(RHO_LC(1,ind)),'*');
if isempty(ind)
    ind = find(p_LC(1,:)<0.005);
    plot(X_TIME(1,ind),RHO_LC(1,ind),'o');
end
xlabel('time from target');

subplot(2,2,4)
plot(X_TIME(2,:),abs(RHO_LC(2,:)))
box off;
hold on;
ylim([0 1])
xlim([-600 800])
ind = find(p_LC(2,:)<0.0001);
plot(X_TIME(2,ind),abs(RHO_LC(2,ind)),'*');
if isempty(ind)
    ind = find(p_LC(2,:)<0.005);
    plot(X_TIME(2,ind),RHO_LC(2,ind),'o');
end
xlabel('time from movement');

suptitle(strcat(NOME(1:8),'-',NOME(10:12),'-RHO-EPOCHS'));

cd(Results_dir)
filename = strcat(NOME,'_RHO_EPOCHS');
print(FFFFFF, '-dpdf', filename, '-r600')




















[a,b] = nanmin(nanmin(P_CN_T_smooth));
[c,d] = nanmin(nanmin(P_CN_M_smooth));

if a<c
    SUGG_EPOCH = [time_T(b)-100;time_T(b)+100];
    SUGG_ALIGN = 'T';
elseif a>c
    SUGG_EPOCH = [time_M(d)-100;time_M(d)+100];
    SUGG_ALIGN = 'M';
end


disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp(strcat('Suggested epoch is ',num2str(SUGG_EPOCH(1)),':',num2str(SUGG_EPOCH(2)),'ms @',SUGG_ALIGN ));
disp('!!! enter the epoch to analyse, here !!!');
disp('!!! START : END !!!');
epoch =(input('Enter here: ','s'));
[C,matches] = strsplit(epoch,{':'});

if ~strcmp(epoch,'0:0')
    EPOCH_START = str2num(cell2mat(C(1)));
    EPOCH_END = str2num(cell2mat(C(2)));
    disp('!!! Data accepted successfully !!!');
    disp('!!! What should I align the spikes to? !!!');
    ALIGN = upper(input('Enter here: ','s'));
    
    
    
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
    
    
    suptitle(strcat(NOME(1:8),'-',NOME(10:12),'-PSTH-SS-PREV-EPOCH-DYN'));
    
    
    cd(Results_dir)
    filename = strcat(NOME,'_PSTH_SS_PREV_EPOCH_DYN');
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
    
    
    
    
    
    %%
    
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
    save(MERGE_file,'ALIGN','XVAL_CW','VALUE_C','RHO_LC_PSTH','p_LC_PSTH','VALUE_C_BEF','VALUE_C_DUR','VALUE_C_AFT','VALUE_C_TOT','per_corr_TOT','XVAL_TOT',...
        'GROUP','EPOCH_TYPE','EPOCH_START','EPOCH_END','-append');
    save(ALLCELLS_file,'ALIGN','XVAL_CW','VALUE_C','RHO_LC_PSTH','p_LC_PSTH','VALUE_C_BEF','VALUE_C_DUR','VALUE_C_AFT','VALUE_C_TOT','per_corr_TOT','XVAL_TOT',...
        'GROUP','EPOCH_TYPE','EPOCH_START','EPOCH_END','-append');
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %
    
    
else
    disp('>>> epoch analysis aborted <<<');
end







%%%%%%%%%%%%%%%%%%%%%%%%% GET THE REAL DELTAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


STRING{1} = 'BASELINE'; STRING{2} = 'TARGET'; STRING{3} = 'MOVEMENT'; STRING{4} = 'POST MOVEMENT';

for IND =1:4
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp(strcat(['!!! enter the',' ', STRING{IND}, ' epoch here !!!']));
    disp('!!! START : END !!!');
    epoch =(input('Enter here: ','s'));
    [C,matches] = strsplit(epoch,{':'});
    
    if ~strcmp(epoch,'0:0')
        DELTA.START(IND) = str2num(cell2mat(C(1)));
        DELTA.END(IND) = str2num(cell2mat(C(2)));
        disp('!!! Data accepted successfully !!!');
    else
        DELTA.START(IND) = NaN;
        DELTA.END(IND) = NaN;
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FOR TILING DELTAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp(strcat('!!! enter the sign for PREV here !!!'));
sign(1) =(input('Enter here: '));

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp(strcat('!!! enter the sign for CURR here !!!'));
sign(2) =(input('Enter here: '));

DELTA.HELPER = sign;






for IND=1:4
    
    EPOCH_START = DELTA.START(IND);
    EPOCH_END = DELTA.END(IND);
    
    if IND==1 | IND==2 ALIGN = 'T'; end
    if IND==3 | IND==4 ALIGN = 'M'; end
    
    if isnan(EPOCH_START) EPOCH_START=0; end
    if isnan(EPOCH_END) EPOCH_END=0; end
    
    epoch = strcat(num2str(EPOCH_START),':',num2str(EPOCH_END));
    
    if ~strcmp(epoch,'0:0')
        Delta_epochs
        DELTA.EpochFlag(IND)=1;
    else
        DELTA.EpochFlag(IND)=0;
        disp('>>> epoch analysis aborted <<<');
    end
    
end





for IND=1:4
    
    EPOCH_START = DELTA.START(IND);
    EPOCH_END = DELTA.END(IND);
    
    if IND==1 | IND==2 ALIGN = 'T'; end
    if IND==3 | IND==4 ALIGN = 'M'; end
    
    if isnan(EPOCH_START) EPOCH_START=0; end
    if isnan(EPOCH_END) EPOCH_END=0; end
    
    epoch = strcat(num2str(EPOCH_START),':',num2str(EPOCH_END));
    
    if ~strcmp(epoch,'0:0')
        FigureMaker_LC_RT_DELTA
    else
        disp('>>> epoch analysis aborted <<<');
    end
    
end




save(POP_file,'DELTA','-append');
save(MERGE_file,'DELTA','-append');
save(ALLCELLS_file,'DELTA','-append');


