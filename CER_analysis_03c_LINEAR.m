

%% function CER_analysis_03c_LINEAR
% Memory analysis

NUM = 19;

%%
% 
WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;



for i=2:length(Infos)
    
    if (Infos(i,10)==1) && (Infos(i-1,10)==2)     % Wrong -> Correct
        Infos(i,16)=1;
    end
    
    if (Infos(i,10)==2) && (Infos(i-1,10)==1)     % Correct -> Wrong
        Infos(i,16)=2;
    end
    
    if (Infos(i,10)==1) && (Infos(i-1,10)==1)     % Correct -> Correcct
        Infos(i,16)=3;
    end
    
    if (Infos(i,10)==2) && (Infos(i-1,10)==2)     % Wrong -> Wrong
        Infos(i,16)=4;
    end
    
end



% Start_T = -450;
% End_T  =1050;
% StartLIM_T = -400; EndLIM_T = 800;
% 
% Start_M = -750;
% End_M  =750;
% StartLIM_M = -500; EndLIM_M = 700;


%%% NEW  %%% 11/27/19
Start_T = -1550;
End_T  =1550;
StartLIM_T = -1400; EndLIM_T = 700;

Start_M = -1550;
End_M  =1550;
StartLIM_M = -700; EndLIM_M = 1400;

FF = figure;
SIG =30;


% Wrong -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==1)+CHANGE-1;
WC_NUM=length(IND);
subplot(2,2,1)
WC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
subplot(2,2,2)
WC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim; 


% Correct -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==2)+CHANGE-1;
CW_NUM=length(IND);
subplot(2,2,1)
CW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim;
subplot(2,2,2)
CW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Correct -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==3)+CHANGE-1;
CC_NUM=length(IND);
subplot(2,2,1)
CC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
subplot(2,2,2)
CC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;

% Wrong -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==4)+CHANGE-1;
WW_NUM=length(IND);
subplot(2,2,1)
WW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim;
subplot(2,2,2)
WW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;



%MAIN TASK
subplot(2,2,1)
IND = CHANGE-15:CHANGE;
MT_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim;
subplot(2,2,2)
MT_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM10 = ylim; 


subplot(2,2,1)
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
hold on;
plot([0 0],ylim,'-k','linewidth',1);
subplot(2,2,2) 
hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
plot([0 0],ylim,'-k','linewidth',1);





subplot(2,2,3)
box off; axis off;
text(0,0.4,strcat('Wrong->Correct =',num2str(WC_NUM)),'color','m');
text(0,0.6,strcat('Correct->Wrong =',num2str(CW_NUM)),'color','b');
text(0,0.8,strcat('Correct->Correct =',num2str(CC_NUM)),'color','g');
text(0,1.0,strcat('Wrong->Wrong =',num2str(WW_NUM)),'color','r');

% text(0,0,strcat('CHANGE = ',num2str(CHANGE),'; LEARNT = ',num2str(LEARNT)),'color','K');





% % suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));
% % 
% % cd(Results_dir)
% % filename = 'PSTH_Simple spike_CUATRO_TODO';
% % print(FF, '-dpdf', filename, '-r400')





% baseline average -------
WN_T_avg = nanmean([WW_T;WC_T]);
WN_T_std = nanstd([WW_T;WC_T]);
CN_T_avg = nanmean([CW_T;CC_T]);
CN_T_std = nanstd([CW_T;CC_T]);
time_T = Start_T:End_T;

% PM average -------
NW_M_avg = nanmean([WW_M;CW_M]);
NW_M_std = nanstd([WW_M;CW_T]);
NC_M_avg = nanmean([WC_M;CC_M]);
NC_M_std = nanstd([WC_M;CC_M]);
time_M = -End_M:750;


% STATS -------------------------

%baseline
ind = find(-400<=time_T & time_T<=-100);
WW_T_stats(1) = nanmean(WW_T(ind));
CW_T_stats(1) = nanmean(CW_T(ind));
CC_T_stats(1) = nanmean(CC_T(ind));
WC_T_stats(1) = nanmean(WC_T(ind));
WW_T_stats(2) = nanstd(WW_T(ind));
CW_T_stats(2) = nanstd(CW_T(ind));
CC_T_stats(2) = nanstd(CC_T(ind));
WC_T_stats(2) = nanstd(WC_T(ind));

CN_T_stats(1) = nanmean(CN_T_avg(ind));
CN_T_stats(2) = nanstd(CN_T_avg(ind));

WN_T_stats(1) = nanmean(WN_T_avg(ind));
WN_T_stats(2) = nanstd(WN_T_avg(ind));

MT_T_stats(1) = nanmean(MT_T(ind));
MT_T_stats(2) = nanstd(MT_T(ind));


%PM
ind = find(400<=time_M & time_M<=700);
WW_M_stats(1) = nanmean(WW_M(ind));
CW_M_stats(1) = nanmean(CW_M(ind));
CC_M_stats(1) = nanmean(CC_M(ind));
WC_M_stats(1) = nanmean(WC_M(ind));
WW_M_stats(2) = nanstd(WW_M(ind));
CW_M_stats(2) = nanstd(CW_M(ind));
CC_M_stats(2) = nanstd(CC_M(ind));
WC_M_stats(2) = nanstd(WC_M(ind));

NC_M_stats(1) = nanmean(NC_M_avg(ind));
NC_M_stats(2) = nanstd(NC_M_avg(ind));

NW_M_stats(1) = nanmean(NW_M_avg(ind));
NW_M_stats(2) = nanstd(NW_M_avg(ind));

MT_M_stats(1) = nanmean(MT_M(ind));
MT_M_stats(2) = nanstd(MT_M(ind));






subplot(2,2,4)
box off;
hold on;
bar(1,WW_T_stats(1),'r')
bar(2,WC_T_stats(1),'m')
bar(3,MT_T_stats(1),'k')
bar(4,CW_T_stats(1),'b')
bar(5,CC_T_stats(1),'g')

bar(8,WW_M_stats(1),'r')
bar(9,CW_M_stats(1),'b')
bar(10,MT_M_stats(1),'k')
bar(11,CC_M_stats(1),'g')
bar(12,WC_M_stats(1),'m')
errorbar([1 2 3 4 5],[WW_T_stats(1)  WC_T_stats(1) MT_T_stats(1) CW_T_stats(1) CC_T_stats(1)],[WW_T_stats(2) WC_T_stats(2) MT_T_stats(2) CW_T_stats(2) CC_T_stats(2)],'.');
errorbar([8 9 10 11 12],[WW_M_stats(1) CW_M_stats(1) MT_M_stats(1) CC_M_stats(1) WC_M_stats(1)],[WW_M_stats(2) CW_M_stats(2) MT_M_stats(2) CC_M_stats(2) WC_M_stats(2)],'.');
xlim([0 13])

plot([6.5 6.5],ylim,'--k')
YLIM = ylim;
ylim([-5 YLIM(2)])



% NW_P = ttest_n(WW_M_stats(1),CW_M_stats(1),WW_M_stats(2),CW_M_stats(2),WW_NUM,CW_NUM);
% NC_P = ttest_n(WC_M_stats(1),CC_M_stats(1),WC_M_stats(2),CC_M_stats(2),WC_NUM,CC_NUM);
% NW_NUM = CW_NUM+WW_NUM; NC_NUM = CC_NUM+WC_NUM;
% NW_P = ttest_n(NW_M_stats(1),NC_M_stats(1),NW_M_stats(2),NC_M_stats(2),NW_NUM,NC_NUM);




suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = 'PSTH_Simple spike_CUATRO_TODO';
print(FF, '-dpdf', filename, '-r400')







save(MERGE_file,'WW_M_stats', 'CW_M_stats', 'CC_M_stats', 'WC_M_stats','WW_T_stats', 'CW_T_stats', 'CC_T_stats', 'WC_T_stats','MT_T_stats', 'MT_M_stats','-append');
save(ALLCELLS_file,'WW_M_stats', 'CW_M_stats', 'CC_M_stats', 'WC_M_stats','WW_T_stats', 'CW_T_stats', 'CC_T_stats', 'WC_T_stats','WC_T_stats','MT_T_stats', 'MT_M_stats','-append');




% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure for D&D stats ----------
% % 
% % 
% % FF = figure;
% % 
% % subplot(4,4,6)
% % box off;
% % hold on;
% % bar(1,WN_T_stats(1),'facecolor',WRG_COLOUR)
% % bar(2,MT_T_stats(1),'facecolor',OT_COLOUR)
% % bar(3,CN_T_stats(1),'facecolor',COR_COLOUR)
% % 
% % 
% % bar(5,NW_M_stats(1),'facecolor',WRG_COLOUR)
% % bar(6,MT_M_stats(1),'facecolor',OT_COLOUR)
% % bar(7,NC_M_stats(1),'facecolor',COR_COLOUR)
% % 
% % errorbar([1 2 3],[WN_T_stats(1) MT_T_stats(1) CN_T_stats(1)],[WN_T_stats(2) MT_T_stats(2) CN_T_stats(2)],'.k');
% % errorbar([5 6 7],[NW_M_stats(1) MT_M_stats(1) NC_M_stats(1)],[NW_M_stats(2) MT_M_stats(2) NC_M_stats(2)],'.k');
% % xlim([0 8])
% % 
% % plot([4 4],ylim,'--k')
% % YLIM = ylim;
% % ylim([-5 YLIM(2)])
% % 
% % 
% % suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));
% % 
% % cd(Results_dir)
% % filename = 'D&D_Stats';
% % print(FF, '-dpdf', filename, '-r400')
% % 








% % % %% ------------------------------------------------------
% % % 
% % % % baseline average -------
% % % WN_T_avg = nanmean([WW_T;WC_T]);
% % % WN_T_std = nanstd([WW_T;WC_T]);
% % % CN_T_avg = nanmean([CW_T;CC_T]);
% % % CN_T_std = nanstd([CW_T;CC_T]);
% % % time_T = Start_T:End_T;
% % % 
% % % % PM average -------
% % % NW_M_avg = nanmean([WW_M;CW_M]);
% % % NW_M_std = nanstd([WW_M;CW_T]);
% % % NC_M_avg = nanmean([WC_M;CC_M]);
% % % NC_M_std = nanstd([WC_M;CC_M]);
% % % time_M = -End_M:750;
% % % 
% % % 
% % % % STATS -------------------------
% % % 
% % % %baseline
% % % ind = find(-400<=time_T & time_T<=-100);
% % % WW_T_stats(1) = nanmean(WW_T(ind));
% % % CW_T_stats(1) = nanmean(CW_T(ind));
% % % CC_T_stats(1) = nanmean(CC_T(ind));
% % % WC_T_stats(1) = nanmean(WC_T(ind));
% % % WW_T_stats(2) = nanstd(WW_T(ind));
% % % CW_T_stats(2) = nanstd(CW_T(ind));
% % % CC_T_stats(2) = nanstd(CC_T(ind));
% % % WC_T_stats(2) = nanstd(WC_T(ind));
% % % 
% % % CN_T_stats(1) = nanmean(CN_T_avg(ind));
% % % CN_T_stats(2) = nanstd(CN_T_avg(ind));
% % % 
% % % WN_T_stats(1) = nanmean(WN_T_avg(ind));
% % % WN_T_stats(2) = nanstd(WN_T_avg(ind));
% % % 
% % % MT_T_stats(1) = nanmean(MT_T(ind));
% % % MT_T_stats(2) = nanstd(MT_T(ind));
% % % 
% % % 
% % % %PM
% % % ind = find(400<=time_M & time_M<=700);
% % % WW_M_stats(1) = nanmean(WW_M(ind));
% % % CW_M_stats(1) = nanmean(CW_M(ind));
% % % CC_M_stats(1) = nanmean(CC_M(ind));
% % % WC_M_stats(1) = nanmean(WC_M(ind));
% % % WW_M_stats(2) = nanstd(WW_M(ind));
% % % CW_M_stats(2) = nanstd(CW_M(ind));
% % % CC_M_stats(2) = nanstd(CC_M(ind));
% % % WC_M_stats(2) = nanstd(WC_M(ind));
% % % 
% % % NC_M_stats(1) = nanmean(NC_M_avg(ind));
% % % NC_M_stats(2) = nanstd(NC_M_avg(ind));
% % % 
% % % NW_M_stats(1) = nanmean(NW_M_avg(ind));
% % % NW_M_stats(2) = nanstd(NW_M_avg(ind));
% % % 
% % % MT_M_stats(1) = nanmean(MT_M(ind));
% % % MT_M_stats(2) = nanstd(MT_M(ind));
% % % 
% % % 
% % % 












% end
