

%% function CER_analysis_03pp_LINEAR
% Memory analysis

NUM = 14;

INFOS = Infos;
Infos(2:end,10) = Infos(1:end-1,10); % 1+N



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



for i=2:length(INFOS)
    
    if (INFOS(i,10)==1) && (INFOS(i-1,10)==2)     % Wrong -> Correct
        Infos(i,17)=1;
    end
    
    if (INFOS(i,10)==2) && (INFOS(i-1,10)==1)     % Correct -> Wrong
        Infos(i,17)=2;
    end
    
    if (INFOS(i,10)==1) && (INFOS(i-1,10)==1)     % Correct -> Correcct
        Infos(i,17)=3;
    end
    
    if (INFOS(i,10)==2) && (INFOS(i-1,10)==2)     % Wrong -> Wrong
        Infos(i,17)=4;
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

Start_R = 200; End_R = 800;
StartLIM_R = 300; EndLIM_R = 700;

Start_T = -1500; End_T  =500;
StartLIM_T = -1400; EndLIM_T = 400;

Start_M = -500; End_M  =400;
StartLIM_M = -400; EndLIM_M = 300;


FF = figure;
SIG =30;


% Wrong -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,17)==1)+CHANGE-1;
subplot(6,6,1);
WC_R = PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_R,End_R,SIG,'m',1,0); % Aligned to target
hold on;
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim;


% Correct -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,17)==2)+CHANGE-1;
subplot(6,6,1)
CW_R = PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_M,End_M,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Correct -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,17)==3)+CHANGE-1;
subplot(6,6,1)
CC_R = PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_M,End_M,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;

% Wrong -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,17)==4)+CHANGE-1;
subplot(6,6,1)
WW_R = PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_M,End_M,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;

xlim([StartLIM_R EndLIM_R]);








% Wrong -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==1)+CHANGE-1;
WC_NUM=length(IND);
subplot(6,6,2)
WC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM1 = ylim;
subplot(6,6,3)
WC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'m',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM2 = ylim;


% Correct -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==2)+CHANGE-1;
CW_NUM=length(IND);
subplot(6,6,2)
CW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM3 = ylim;
subplot(6,6,3)
CW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'b',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM4 = ylim;


% Correct -> Correct
IND = find(Infos(CHANGE:CHANGE+NUM,16)==3)+CHANGE-1;
CC_NUM=length(IND);
subplot(6,6,2)
CC_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM5 = ylim;
subplot(6,6,3)
CC_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'g',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM6 = ylim;

% Wrong -> Wrong
IND = find(Infos(CHANGE:CHANGE+NUM,16)==4)+CHANGE-1;
WW_NUM=length(IND);
subplot(6,6,2)
WW_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM7 = ylim;
subplot(6,6,3)
WW_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'r',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM8 = ylim;


%MAIN TASK
subplot(6,6,2)
IND = CHANGE-15:CHANGE;
MT_T = PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); xlabel([]); xlim([StartLIM_T EndLIM_T]); YLIM9 = ylim; set(gca,'YTick',[]); ylabel([]); xlabel([]);
subplot(6,6,3)
MT_M = PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); set(gca,'YTick',[]); ylabel([]); xlabel([]); xlim([StartLIM_M EndLIM_M]); YLIM10 = ylim;
subplot(6,6,1)
MT_R = PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_R,End_R,SIG,'k',1,0); % Aligned to target
set(gca,'fontsize',7); xlim([StartLIM_R EndLIM_R]);


subplot(6,6,2)
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
hold on;
plot([0 0],ylim,'-k','linewidth',1);
subplot(6,6,3)
hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
plot([0 0],ylim,'-k','linewidth',1);
subplot(6,6,1)
hold on;
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);



clear temp1 temp2 temp3 temp4 temp5 temp6;
temp1 = NaN(1,2001); temp2 = NaN(1,2001);
temp3 = NaN(1,901); temp4 = NaN(1,901);
temp5 = NaN(1,601); temp6 = NaN(1,601);

clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==4)+CHANGE-1;
try temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,17)==4)+CHANGE-1;
try temp5 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end

clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==2)+CHANGE-1;
try temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,17)==2)+CHANGE-1;
try temp6 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end


WN_T = [temp1;temp2]; WN_M = [temp3;temp4];
% WN = [WN_T WN_M];
% WN_T = [temp1;temp3]; WN_M = [temp2;temp4];
WN_R = [temp5; temp6];
try
    clear WN; WN = [WN_R WN_T WN_M];
catch
    try
        clear WN; WN = [WN_R WN_T(2:end,:) WN_M(2:end,:)];
    catch
        clear WN; WN = [WN_R(2:end,:) WN_T WN_M];
    end
end


clear temp1 temp2 temp3 temp4 temp5 temp6;
temp1 = NaN(1,2001); temp2 = NaN(1,2001);
temp3 = NaN(1,901); temp4 = NaN(1,901);
temp5 = NaN(1,601); temp6 = NaN(1,601);

clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==1)+CHANGE-1;
try temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,17)==1)+CHANGE-1;
try temp5 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end

clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,16)==3)+CHANGE-1;
try temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE:CHANGE+NUM,17)==3)+CHANGE-1;
try temp6 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end


CN_T = [temp1;temp2]; CN_M = [temp3;temp4];
% CN = [CN_T CN_M];
% CN_T = [temp1;temp3]; CN_M = [temp2;temp4];
CN_R = [temp5; temp6];
try
    clear CN; CN = [CN_R CN_T CN_M];
catch
    try
        clear CN; CN = [CN_R(2:end,:) CN_T CN_M];
    catch
        clear CN; CN = [CN_R CN_T(2:end,:) CN_M(2:end,:)];
    end
end

% CN_T = [temp1;temp2]; CN_M = [temp3;temp4]; CN = [CN_T CN_M];












FF = figure;

IND = CHANGE-12:CHANGE;

subplot(6,29,[1 4]); hold on;
errorline_n(Start_R:End_R,nanmean(WN_R),nanstd(WN_R)/sqrt(size(WN_R,1)),1,[1 0 0]);
errorline_n(Start_R:End_R,nanmean(CN_R),nanstd(CN_R)/sqrt(size(CN_R,1)),1,[0 0 1]);
PSTH_n(Spikes.S(IND,:),INFOS(IND,11),Start_R,End_R,SIG,'k',1,0); 
xlim([StartLIM_R EndLIM_R]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
% ylim([0 80]);

subplot(6,29,[5 22]); hold on;
errorline_n(Start_T:End_T,nanmean(WN_T),nanstd(WN_T)/sqrt(size(WN_T,1)),1,[1 0 0]);
errorline_n(Start_T:End_T,nanmean(CN_T),nanstd(CN_T)/sqrt(size(CN_T,1)),1,[0 0 1]);
PSTH_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG,'k',1,0);
xlim([StartLIM_T EndLIM_T]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]); set(gca,'YTick',[]); ylabel([]);

subplot(6,29,[23 29]); hold on;
errorline_n(Start_M:End_M,nanmean(WN_M),nanstd(WN_M)/sqrt(size(WN_M,1)),1,[1 0 0]);
errorline_n(Start_M:End_M,nanmean(CN_M),nanstd(CN_M)/sqrt(size(CN_M,1)),1,[0 0 1]);
PSTH_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,'k',1,0);
xlim([StartLIM_M EndLIM_M]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]); set(gca,'YTick',[]); ylabel([]);





alltime = [Start_R: End_R Start_T:End_T Start_M:End_M];
LIM1 = find(alltime==StartLIM_R,1);
LIM2 = find(alltime==EndLIM_R,1);
LIM3 = find(alltime==StartLIM_T,1);
LIM4 = find(alltime==EndLIM_T,2); LIM4=LIM4(2);
LIM5 = find(alltime==StartLIM_M); LIM5=LIM5(2);
LIM6 = find(alltime==EndLIM_M); LIM6 = LIM6(3);


clear P_ttest
parfor ii=1:size(WN,2)
    P_ttest(1,ii) = ttest_NN(WN(:,ii),CN(:,ii));
end

S1 = subplot(6,29,[1 4]+29); hold on; hold on;
plot(StartLIM_R:EndLIM_R,log10(P_ttest(LIM1:LIM2)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_R EndLIM_R]);
S2 = subplot(6,29,[5 22]+29); hold on;
plot(StartLIM_T:EndLIM_T,log10(P_ttest(LIM3:LIM4)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_T EndLIM_T]); set(gca,'YTick',[]); ylabel([]);
S3 = subplot(6,29,[23 29]+29); hold on;
plot(StartLIM_M:EndLIM_M,log10(P_ttest(LIM5:LIM6)));
plot(xlim,[log10(0.05) log10(0.05)],':r');
xlim([StartLIM_M EndLIM_M]); set(gca,'YTick',[]); ylabel([]);



%%%%%%%%%%%%%%%%%%% NOW ANOVA %%%%%%%%%%%%%%%%%%



% % % % % IND1 = CHANGE-20:CHANGE;
% % % % % IND2 = CHANGE:CHANGE+20;
% % % % % IND3 = CHANGE+21:CHANGE+40;
% % % % % IND4 = CHANGE+40:nanmin([LEARNT CHANGE+60 size(Infos,1)]);
% % % % % if isempty(IND4)
% % % % %     clear IND3 IND4;
% % % % %     IND4 = nanmin([LEARNT CHANGE+60 size(Infos,1)])-10:nanmin([LEARNT CHANGE+60 size(Infos,1)]);
% % % % %     IND3 = nanmin([LEARNT CHANGE+60 size(Infos,1)])-20:nanmin([LEARNT CHANGE+60 size(Infos,1)])-10;
% % % % % end
% % % % %
% % % % % IND1 = IND1(find(Infos(IND1,10)==1));
% % % % % IND2 = IND2(find(Infos(IND2,10)==1));
% % % % % IND3 = IND3(find(Infos(IND3,10)==1));
% % % % % IND4 = IND4(find(Infos(IND4,10)==1));
% % % % %
% % % % %
% % % % %
% % % % %
% % % % %
% % % % %
% % % % % subplot(6,29,[1 4]+29*3); hold on;
% % % % % PSTH_n(Spikes.S(IND1,:),Infos(IND1,11),Start_R,End_R,SIG,'k',1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND2,:),Infos(IND2,11),Start_R,End_R,SIG,[0.3 0.3 0.3],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND3,:),Infos(IND3,11),Start_R,End_R,SIG,[0.6 0.6 0.6],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND4,:),Infos(IND4,11),Start_R,End_R,SIG,[0.9 0.9 0.9],1,0); hold on;
% % % % % ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
% % % % % plot([0 0],ylim,'-k','linewidth',1); xlim([StartLIM_R EndLIM_R]);
% % % % %
% % % % %
% % % % % subplot(6,29,[5 22]+29*3); hold on;
% % % % % PSTH_n(Spikes.S(IND1+1,:),Infos(IND1+1,4),Start_T,End_T,SIG,'k',1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND2+1,:),Infos(IND2+1,4),Start_T,End_T,SIG,[0.3 0.3 0.3],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND3+1,:),Infos(IND3+1,4),Start_T,End_T,SIG,[0.6 0.6 0.6],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND4+1,:),Infos(IND4+1,4),Start_T,End_T,SIG,[0.9 0.9 0.9],1,0); hold on;
% % % % % ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
% % % % % plot([0 0],ylim,'-k','linewidth',1); xlim([StartLIM_T EndLIM_T]); set(gca,'YTick',[]); ylabel([]);
% % % % %
% % % % %
% % % % % subplot(6,29,[23 29]+29*3); hold on;
% % % % % PSTH_n(Spikes.S(IND1+1,:),Infos(IND1+1,11),Start_M,End_M,SIG,'k',1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND2+1,:),Infos(IND2+1,11),Start_M,End_M,SIG,[0.3 0.3 0.3],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND3+1,:),Infos(IND3+1,11),Start_M,End_M,SIG,[0.6 0.6 0.6],1,0); % Aligned to target
% % % % % PSTH_n(Spikes.S(IND4+1,:),Infos(IND4+1,11),Start_M,End_M,SIG,[0.9 0.9 0.9],1,0); hold on;
% % % % % ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
% % % % % plot([0 0],ylim,'-k','linewidth',1); xlim([StartLIM_M EndLIM_M]); set(gca,'YTick',[]); ylabel([]); xlabel([]);
% % % % %
% % % % %
% % % % %
% % % % %
% % % % % P1_R = PSTH_RETURN_n(Spikes.S(IND1,:),Infos(IND1,4),Start_R,End_R,SIG); % Aligned to target
% % % % % P2_R = PSTH_RETURN_n(Spikes.S(IND2,:),Infos(IND2,4),Start_R,End_R,SIG); % Aligned to target
% % % % % P3_R = PSTH_RETURN_n(Spikes.S(IND3,:),Infos(IND3,4),Start_R,End_R,SIG); % Aligned to target
% % % % % P4_R = PSTH_RETURN_n(Spikes.S(IND4,:),Infos(IND4,4),Start_R,End_R,SIG); % Aligned to target
% % % % %
% % % % % P1_T = PSTH_RETURN_n(Spikes.S(IND1+1,:),Infos(IND1+1,4),Start_T,End_T,SIG); % Aligned to target
% % % % % P2_T = PSTH_RETURN_n(Spikes.S(IND2+1,:),Infos(IND2+1,4),Start_T,End_T,SIG); % Aligned to target
% % % % % P3_T = PSTH_RETURN_n(Spikes.S(IND3+1,:),Infos(IND3+1,4),Start_T,End_T,SIG); % Aligned to target
% % % % % P4_T = PSTH_RETURN_n(Spikes.S(IND4+1,:),Infos(IND4+1,4),Start_T,End_T,SIG); % Aligned to target
% % % % %
% % % % % P1_M = PSTH_RETURN_n(Spikes.S(IND1+1,:),Infos(IND1+1,11),Start_M,End_M,SIG); % Aligned to target
% % % % % P2_M = PSTH_RETURN_n(Spikes.S(IND2+1,:),Infos(IND2+1,11),Start_M,End_M,SIG); % Aligned to target
% % % % % P3_M = PSTH_RETURN_n(Spikes.S(IND3+1,:),Infos(IND3+1,11),Start_M,End_M,SIG); % Aligned to target
% % % % % P4_M = PSTH_RETURN_n(Spikes.S(IND4+1,:),Infos(IND4+1,11),Start_M,End_M,SIG); % Aligned to target
% % % % %
% % % % %
% % % % % P1 = [P1_R P1_T P1_M]; P2 = [P2_R P2_T P2_M]; P3 = [P3_R P3_T P3_M]; P4 = [P4_R P4_T P4_M];
% % % % % LENN = nanmax([size(P1,1) size(P2,1) size(P3,1) size(P4,1)]);


% % % % % parfor ii=1:size(P1,2)
% % % % %     ANOVA_MAT = NaN(LENN,4);
% % % % %     ANOVA_MAT(1:length(P1(:,ii)),1) = P1(:,ii);
% % % % %     ANOVA_MAT(1:length(P2(:,ii)),2) = P2(:,ii);
% % % % %     ANOVA_MAT(1:length(P3(:,ii)),3) = P3(:,ii);
% % % % %     ANOVA_MAT(1:length(P4(:,ii)),4) = P4(:,ii);
% % % % %     P_anova(1,ii) = anova1(ANOVA_MAT);
% % % % %     delete(gcf); delete(gcf);
% % % % % end
% % % % %
% % % % %
% % % % %
% % % % %
% % % % % S4 = subplot(6,29,[1 4]+29*4); hold on;
% % % % % plot(StartLIM_R:EndLIM_R,log10(P_anova(LIM1:LIM2)));
% % % % % plot(xlim,[log10(0.05) log10(0.05)],':r');
% % % % % xlim([StartLIM_R EndLIM_R]);
% % % % %
% % % % % S5 = subplot(6,29,[5 22]+29*4); hold on;
% % % % % plot(StartLIM_T:EndLIM_T,log10(P_anova(LIM3:LIM4)));
% % % % % plot(xlim,[log10(0.05) log10(0.05)],':r');
% % % % % xlim([StartLIM_T EndLIM_T]); set(gca,'YTick',[]); ylabel([]);
% % % % %
% % % % % S6 = subplot(6,29,[23 29]+29*4); hold on;
% % % % % plot(StartLIM_M:EndLIM_M,log10(P_anova(LIM5:LIM6)));
% % % % % plot(xlim,[log10(0.05) log10(0.05)],':r');
% % % % % xlim([StartLIM_M EndLIM_M]); set(gca,'YTick',[]); ylabel([]);
% % % linkaxes([S1 S2 S3 S4 S5 S6],'y');


linkaxes([S1 S2 S3],'y');








% % P_anova(isnan(P_anova))=1;
% % P_ttest(isnan(P_ttest))=1;
% % CORR_VAL = corr(P_anova(:),P_ttest(:));

% % STATE_DELTA.P_anova = P_anova;
STATE_DELTA.P_ttest = P_ttest;
% % STATE_DELTA.CORR = CORR_VAL;

suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));

cd(Results_dir)
filename = 'PSTH_Simple spike_CUATRO_STATECHANGE';
print(FF, '-dpdf', filename, '-r400')

cd('E:\NAVEEN_Work\Cerebellum\Data\deltas');
filename = strcat(NOME,'_delta');
print(FF, '-dpdf', filename, '-r400');


save(MERGE_file,'STATE_DELTA','-append');
% save(ALLCELLS_file,'STATE_DELTA','-append');
save(POP_file,'STATE_DELTA','-append');







%% %%%%%%%%%% DELTA FOR OT 

clear temp1 temp2 temp3 temp4 temp5 temp6;
temp1 = NaN(1,2001); temp2 = NaN(1,2001);
temp3 = NaN(1,901); temp4 = NaN(1,901);
temp5 = NaN(1,601); temp6 = NaN(1,601);

NUM = CHANGE-1;

clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,16)==4)+CHANGE-1-NUM;
try temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,17)==4)+CHANGE-1-NUM;
try temp5 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end

clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,16)==2)+CHANGE-1-NUM;
try temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,17)==2)+CHANGE-1-NUM;
try temp6 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end


WN_T = [temp1;temp2]; WN_M = [temp3;temp4];
% WN = [WN_T WN_M];
% WN_T = [temp1;temp3]; WN_M = [temp2;temp4];
WN_R = [temp5; temp6];
try
    clear WN; WN = [WN_R WN_T WN_M];
catch
    try
        clear WN; WN = [WN_R WN_T(2:end,:) WN_M(2:end,:)];
    catch
        clear WN; WN = [WN_R(2:end,:) WN_T WN_M];
    end
end


clear temp1 temp2 temp3 temp4 temp5 temp6;
temp1 = NaN(1,2001); temp2 = NaN(1,2001);
temp3 = NaN(1,901); temp4 = NaN(1,901);
temp5 = NaN(1,601); temp6 = NaN(1,601);

clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,16)==1)+CHANGE-1-NUM;
try temp1 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp3 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,17)==1)+CHANGE-1-NUM;
try temp5 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end

clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,16)==3)+CHANGE-1-NUM;
try temp2 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,4),Start_T,End_T,SIG); catch end
try temp4 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG); catch end
clear IND; IND = find(Infos(CHANGE-NUM:CHANGE,17)==3)+CHANGE-1-NUM;
try temp6 = PSTH_RETURN_n(Spikes.S(IND,:),Infos(IND,11),Start_R,End_R,SIG); catch end


CN_T = [temp1;temp2]; CN_M = [temp3;temp4];
% CN = [CN_T CN_M];
% CN_T = [temp1;temp3]; CN_M = [temp2;temp4];
CN_R = [temp5; temp6];
try
    clear CN; CN = [CN_R CN_T CN_M];
catch
    try
        clear CN; CN = [CN_R(2:end,:) CN_T CN_M];
    catch
        clear CN; CN = [CN_R CN_T(2:end,:) CN_M(2:end,:)];
    end
end

% CN_T = [temp1;temp2]; CN_M = [temp3;temp4]; CN = [CN_T CN_M];




FF = figure;


subplot(6,29,[1 4]); hold on;
errorline_n(Start_R:End_R,nanmean(WN_R),nanstd(WN_R)/sqrt(size(WN_R,1)),1,[1 0 0]);
errorline_n(Start_R:End_R,nanmean(CN_R),nanstd(CN_R)/sqrt(size(CN_R,1)),1,[0 0 1]);
xlim([StartLIM_R EndLIM_R]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]);
% ylim([0 80]);

subplot(6,29,[5 22]); hold on;
errorline_n(Start_T:End_T,nanmean(WN_T),nanstd(WN_T)/sqrt(size(WN_T,1)),1,[1 0 0]);
errorline_n(Start_T:End_T,nanmean(CN_T),nanstd(CN_T)/sqrt(size(CN_T,1)),1,[0 0 1]);
xlim([StartLIM_T EndLIM_T]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]); set(gca,'YTick',[]); ylabel([]);

subplot(6,29,[23 29]); hold on;
errorline_n(Start_M:End_M,nanmean(WN_M),nanstd(WN_M)/sqrt(size(WN_M,1)),1,[1 0 0]);
errorline_n(Start_M:End_M,nanmean(CN_M),nanstd(CN_M)/sqrt(size(CN_M,1)),1,[0 0 1]);
xlim([StartLIM_M EndLIM_M]);
ylim([0 nanmax([YLIM1(2);YLIM2(2);YLIM3(2);YLIM4(2);YLIM5(2);YLIM6(2);YLIM7(2);YLIM8(2);YLIM9(2);YLIM10(2)])]); set(gca,'YTick',[]); ylabel([]);





alltime = [Start_R: End_R Start_T:End_T Start_M:End_M];
LIM1 = find(alltime==StartLIM_R,1);
LIM2 = find(alltime==EndLIM_R,1);
LIM3 = find(alltime==StartLIM_T,1);
LIM4 = find(alltime==EndLIM_T,2); LIM4=LIM4(2);
LIM5 = find(alltime==StartLIM_M); LIM5=LIM5(2);
LIM6 = find(alltime==EndLIM_M); LIM6 = LIM6(3);

% 
% clear P_ttest
% parfor ii=1:size(WN,2)
%     P_ttest(1,ii) = ttest_NN(WN(:,ii),CN(:,ii));
% end
% 
% S1 = subplot(6,29,[1 4]+29); hold on; hold on;
% plot(StartLIM_R:EndLIM_R,log10(P_ttest(LIM1:LIM2)));
% plot(xlim,[log10(0.05) log10(0.05)],':r');
% xlim([StartLIM_R EndLIM_R]);
% S2 = subplot(6,29,[5 22]+29); hold on;
% plot(StartLIM_T:EndLIM_T,log10(P_ttest(LIM3:LIM4)));
% plot(xlim,[log10(0.05) log10(0.05)],':r');
% xlim([StartLIM_T EndLIM_T]); set(gca,'YTick',[]); ylabel([]);
% S3 = subplot(6,29,[23 29]+29); hold on;
% plot(StartLIM_M:EndLIM_M,log10(P_ttest(LIM5:LIM6)));
% plot(xlim,[log10(0.05) log10(0.05)],':r');
% xlim([StartLIM_M EndLIM_M]); set(gca,'YTick',[]); ylabel([]);
% 



linkaxes([S1 S2 S3],'y');







suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12)));



cd('E:\NAVEEN_Work\Cerebellum\Data\deltas');
filename = strcat(NOME,'_delta_OT');
print(FF, '-dpdf', filename, '-r400');
















Infos = INFOS;

% end
