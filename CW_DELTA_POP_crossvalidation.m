% CW_DELTA_POP_crossvalidation




NUM = 29;

WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;


ISWRNG_PREV = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==2)
        ISWRNG_PREV(i,1)=1;
    end
end

ISCORR_PREV = NaN(length(Infos),1);
for i=2:length(Infos)
    if (Infos(i-1,10)==1)
        ISCORR_PREV(i,1)=1;
    end
end

ISWRNG_CURR = NaN(length(Infos),1);
for i=1:length(Infos)
    if (Infos(i,10)==2)
        ISWRNG_CURR(i,1)=1;
    end
end

ISCORR_CURR = NaN(length(Infos),1);
for i=1:length(Infos)
    if (Infos(i,10)==1)
        ISCORR_CURR(i,1)=1;
    end
end


% %%%% OLD
% Start_T = -450;
% End_T  =1050;
% StartLIM_T = -400; EndLIM_T = 800;
% 
% Start_M = -750;
% End_M  =750;
% StartLIM_M = -500; EndLIM_M = 700;

% %%% NEW  %%% 9/9/19
% Start_T = -950;
% End_T  =950;
% StartLIM_T = -800; EndLIM_T = 600;
% 
% Start_M = -950;
% End_M  =950;
% StartLIM_M = -600; EndLIM_M = 700;

%%% NEW  %%% 11/27/19
Start_T = -1550;
End_T  =1550;
StartLIM_T = -1400; EndLIM_T = 700;

Start_M = -1550;
End_M  =1550;
StartLIM_M = -700; EndLIM_M = 1400;



FF = figure;
SIG =40;


% Prev Wrng
IND = find(ISWRNG_PREV(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
IND = IND(randperm(length(IND)));
IND = IND(1:round(length(IND)/2));
W_PREV_NUM=length(IND);
subplot(2,2,1)
W_PREV_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
subplot(2,2,2)
W_PREV_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim; 


% Prev Corr
IND = find(ISCORR_PREV(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
IND = IND(randperm(length(IND)));
IND = IND(1:round(length(IND)/2));
C_PREV_NUM=length(IND);
subplot(2,2,1)
C_PREV_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim; ylabel('Prev','fontsize',12);
subplot(2,2,2)
C_PREV_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Curr Wrng
IND = find(ISWRNG_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
IND = IND(randperm(length(IND)));
IND = IND(1:round(length(IND)/2));
W_CURR_NUM=length(IND);
subplot(2,2,3)
W_CURR_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
subplot(2,2,4)
W_CURR_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]);ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim; 


% Curr Corr
IND = find(ISCORR_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
IND = IND(randperm(length(IND)));
IND = IND(1:round(length(IND)/2));
C_CURR_NUM=length(IND);
subplot(2,2,3)
C_CURR_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim; ylabel('Curr','fontsize',12);
subplot(2,2,4)
C_CURR_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;



% 
% %MAIN TASK
% subplot(2,2,1)
% IND = CHANGE-15:CHANGE;
% MT_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
% set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Prev','fontsize',12);
% subplot(2,2,2)
% MT_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
% set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim; 
% 
% subplot(2,2,3)
% IND = CHANGE-15:CHANGE;
% MT_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
% set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Curr','fontsize',12);
% subplot(2,2,4)
% MT_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
% set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim; 
% 
% 
% 
% 
% for i=1:4
%     subplot(2,2,i)
%     ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
%     hold on;
%     plot([0 0],ylim,'-k','linewidth',1);
% end
% 
% 
% 
% 
% % 
% % 
% % suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));
% % 
% % cd(Results_dir)
% % filename = 'PSTH_Simple spike_DOSyDOS_TODO';
% % print(FF, '-dpdf', filename, '-r400')
% 
% 





%% DELTA FIGURE PLUGIN ---------------------

% % clear CW_DELTA

CW_DELTA.TIME = [Start_T:End_T Start_M:End_M];
PREV_TEMP = (smooth([W_PREV_T(1,:)-C_PREV_T(1,:) W_PREV_M(1,:)-C_PREV_M(1,:)],0.45,'rloess'));
CURR_TEMP = (smooth([W_CURR_T(1,:)-C_CURR_T(1,:) W_CURR_M(1,:)-C_CURR_M(1,:)],0.45,'rloess'));
TEMPTEMP = norm_pos2neg_n([PREV_TEMP; CURR_TEMP]);
CW_DELTA.PREV2 = TEMPTEMP(1:length(PREV_TEMP));
CW_DELTA.CURR2 = TEMPTEMP(length(PREV_TEMP)+1:length(TEMPTEMP));


save(POP_file,'CW_DELTA','-append');









% PREV_TEMP = ([W_PREV_T(1,:)-C_PREV_T(1,:) W_PREV_M(1,:)-C_PREV_M(1,:)]);
% CURR_TEMP = ([W_CURR_T(1,:)-C_CURR_T(1,:) W_CURR_M(1,:)-C_CURR_M(1,:)]);
% TEMPTEMP = norm_pos2neg_n([PREV_TEMP CURR_TEMP]);
% CW_DELTA_PREV = TEMPTEMP(1:length(PREV_TEMP));
% CW_DELTA_CURR = TEMPTEMP(length(PREV_TEMP)+1:length(TEMPTEMP));


% F = figure;
% subplot(10,1,[1 2 3 4])
% hold on;
% plot(1:length(CW_DELTA.TIME),CW_DELTA_PREV,'color',[0.6 0.6 0.6])
% plot(1:length(CW_DELTA.TIME),CW_DELTA.PREV,'linewidth',2)
% value = [1,abs(Start_T)+1,End_T-Start_T-100,End_T-Start_T+100, End_T-Start_T+abs(Start_M), length(CW_DELTA.TIME) ];
% set(gca, 'XTick', value);
% set(gca, 'XTickLabel', {num2str(Start_T),'0',num2str(End_T-100),num2str(Start_M+100),'0',num2str(End_M)},'fontsize',5);
% xlim([1 length(CW_DELTA.TIME)])
% 
% subplot(10,1,[5 6 7 8])
% hold on;
% plot(1:length(CW_DELTA.TIME),CW_DELTA_CURR,'color',[0.6 0.6 0.6])
% plot(1:length(CW_DELTA.TIME),CW_DELTA.CURR,'linewidth',2)
% value = [1,abs(Start_T)+1,End_T-Start_T-100,End_T-Start_T+100, End_T-Start_T+abs(Start_M), length(CW_DELTA.TIME) ];
% set(gca, 'XTick', value);
% set(gca, 'XTickLabel', {num2str(Start_T),'0',num2str(End_T-100),num2str(Start_M+100),'0',num2str(End_M)},'fontsize',5);
% xlim([1 length(CW_DELTA.TIME)])
% 
% 
% subplot(10,1,[9])
% image(CW_DELTA.PREV','CDataMapping','scaled')
% value = [1,abs(Start_T)+1,End_T-Start_T-100,End_T-Start_T+100, End_T-Start_T+abs(Start_M), length(CW_DELTA.TIME) ];
% set(gca, 'XTick', value);
% set(gca, 'XTickLabel', {num2str(Start_T),'0',num2str(End_T-100),num2str(Start_M+100),'0',num2str(End_M)},'fontsize',5);
% xlim([1 length(CW_DELTA.TIME)])
% ax = gca;
% ax.CLim = [-1 1];
% colormap(jet)
% 
% subplot(10,1,[10])
% image(CW_DELTA.CURR','CDataMapping','scaled')
% value = [1,abs(Start_T)+1,End_T-Start_T-100,End_T-Start_T+100, End_T-Start_T+abs(Start_M), length(CW_DELTA.TIME) ];
% set(gca, 'XTick', value);
% set(gca, 'XTickLabel', {num2str(Start_T),'0',num2str(End_T-100),num2str(Start_M+100),'0',num2str(End_M)},'fontsize',5);
% xlim([1 length(CW_DELTA.TIME)])
% ax = gca;
% ax.CLim = [-1 1];
% colormap(jet)
% 
% suptitle(NOME)
% 
% cd(Results_dir)
% filename = strcat(NOME,'ColorDELTA');
% print(F, '-dpdf', filename, '-r400')


disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUCCESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
