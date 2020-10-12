


% function CER_analysis_03d_LINEAR
% dos y dos analysis for OT only




NUM = 19;

WRG_COLOUR = [255 128 128]/255;
COR_COLOUR = [0 179 179]/255;
OT_COLOUR = [81 81 81]/255;


ISWRNG_PREV_OT = NaN(length(Infos),1);
for i=2:CHANGE
    if (Infos(i-1,10)==2)
        ISWRNG_PREV_OT(i,1)=1;
    end
end

ISCORR_PREV_OT = NaN(length(Infos),1);
for i=2:CHANGE
    if (Infos(i-1,10)==1)
        ISCORR_PREV_OT(i,1)=1;
    end
end

ISWRNG_CURR_OT = NaN(length(Infos),1);
for i=1:CHANGE
    if (Infos(i,10)==2)
        ISWRNG_CURR_OT(i,1)=1;
    end
end

ISCORR_CURR_OT = NaN(length(Infos),1);
for i=1:CHANGE
    if (Infos(i,10)==1)
        ISCORR_CURR_OT(i,1)=1;
    end
end


Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;




FF = figure;
SIG =40;


% Prev Wrng
IND = find(ISWRNG_PREV_OT==1);
W_PREV_NUM=length(IND);
subplot(2,2,1)
W_PREV_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
subplot(2,2,2)
W_PREV_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim; 


% Prev Corr
IND = find(ISCORR_PREV_OT==1);
C_PREV_NUM=length(IND);
subplot(2,2,1)
C_PREV_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim; ylabel('Prev','fontsize',12);
subplot(2,2,2)
C_PREV_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Curr Wrng
IND = find(ISWRNG_CURR_OT==1);
W_CURR_NUM=length(IND);
subplot(2,2,3)
W_CURR_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
subplot(2,2,4)
W_CURR_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]);ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim; 


% Curr Corr
IND = find(ISCORR_CURR_OT==1);
C_CURR_NUM=length(IND);
subplot(2,2,3)
C_CURR_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim; ylabel('Curr','fontsize',12);
subplot(2,2,4)
C_CURR_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;




%MAIN TASK
subplot(2,2,1)
IND = 1:CHANGE;
MT_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Prev','fontsize',12);
subplot(2,2,2)
MT_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim; 

subplot(2,2,3)
IND = 1:CHANGE;
MT_T = PSTHe_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,OT_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; ylabel('Curr','fontsize',12);
subplot(2,2,4)
MT_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,OT_COLOUR,1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM10 = ylim; 




for i=1:4
    subplot(2,2,i)
    ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
    hold on;
    plot([0 0],ylim,'-k','linewidth',1);
end






suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-OT'));

cd(Results_dir)
filename = 'PSTH_Simple spike_DOSyDOS_TODO_OT';
print(FF, '-dpdf', filename, '-r400')







clear C_CURR_M W_CURR_M

FF = figure;
SIG =60;

subplot(2,2,4)
hold on;
IND = find(ISCORR_CURR_OT==1);
C_CURR_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,[0 0 1],0.6,0); 
clear IND
IND = find(ISWRNG_CURR_OT==1);
W_CURR_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,[1 0 0],0.6,0); 

% errorline_n(1:size(C_CURR_M,2),C_CURR_M(1,:),C_CURR_M(3,:),2,[0 0 1],0.3,0,0.6)
% errorline_n(1:size(W_CURR_M,2),W_CURR_M(1,:),W_CURR_M(3,:),2,[1 0 0],0.3,0,0.6)
% 
% xlim([1 size(C_CURR_M,2)]);
% ylim([-20 90]);

hold on;
plot([Start_M  Start_M+200],[0 0],'-k')
plot([Start_M Start_M],[0 20],'-k')
axis off;

plot(0,0,'ok','markersize',4);

cd(Results_dir)
filename = 'DOSyDOS_TODO_OT_figure';
print(FF, '-dpdf', filename, '-r400')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating activity in delta epoch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:4
    CW_DELTA.OT_W(index) = NaN;
    CW_DELTA.OT_C(index) = NaN;
end

for index=1:2
    if DELTA.EpochFlag(index)==1
        Time = Start_T:End_T;
        clear INDS
        INDS = find(DELTA.START(index)+95<=Time & Time<=DELTA.END(index)-95);
        % avg ----------------------------------------------
        CW_DELTA.OT_W(index) = nanmean(W_PREV_T(1,INDS));
        CW_DELTA.OT_C(index) = nanmean(C_PREV_T(1,INDS)); 
    end
end

index=3;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+0<=Time & Time<=DELTA.END(index)-0);
    % avg ----------------------------------------------
    CW_DELTA.OT_W(index) = nanmean(W_PREV_M(1,INDS));
    CW_DELTA.OT_C(index) = nanmean(C_PREV_M(1,INDS));
end


index=4;
if DELTA.EpochFlag(index)==1
    Time = Start_M:End_M;
    clear INDS
    INDS = find(DELTA.START(index)+0<=Time & Time<=DELTA.END(index)-0);
    % avg ------------------------------------------
    CW_DELTA.OT_W(index) = nanmean(W_CURR_M(1,INDS));
    CW_DELTA.OT_C(index) = nanmean(C_CURR_M(1,INDS)); 
end



save(POP_file,'CW_DELTA','-append');
save(MERGE_file,'CW_DELTA','-append');
save(ALLCELLS_file,'CW_DELTA','-append');








%%%%% figure 


clear C_CURR_M W_CURR_M

FF = figure;
SIG =60;

subplot(2,2,4)
hold on;

IND = find(ISCORR_CURR_OT==1);
C_CURR_NUM=length(IND);
for i=1:length(IND)
    C_CURR_M(i,:) = PSTH_ONE_n(Spikes.S{IND(i),:},Infos(IND(i),11),Start_M,End_M,SIG,[0 0 1],0.6,0); % Aligned to target
end

IND = find(ISWRNG_CURR_OT==1);
W_CURR_NUM=length(IND);

for i=1:length(IND)
    W_CURR_M(i,:) = PSTH_ONE_n(Spikes.S{IND(i),:},Infos(IND(i),11),Start_M,End_M,SIG,[1 0 0],0.6,0); % Aligned to target
end

hold on;
plot([Start_M  Start_M+200],[0 0],'-k')
plot([Start_M Start_M],[0 20],'-k')
axis off;

plot(0,0,'ok','markersize',4);

cd(Results_dir)
filename = 'DOSyDOS_TODO_OT_figure2';
print(FF, '-dpdf', filename, '-r400')






%%%%%% figure 





