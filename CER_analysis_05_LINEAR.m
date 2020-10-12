%% Part of a batch process
%% Analyses of latency of reward

% Created by NAVEEN ON 10/4/17 at CUMC

% function CER_analysis_05_LINEAR


if isnan(CHANGE)
    CHANGE=1;
end



NUM = 19;
SIG  =30;

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



Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;


% Curr Wrng
IND_NW = find(ISWRNG_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
W_CURR_NUM=length(IND_NW);
W_CURR_M = PSTH_RETURN_n(Spikes.S(IND_NW,:),Infos(IND_NW,11),Start_M,End_M,SIG);

% Curr Corr
IND_NC = find(ISCORR_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
C_CURR_NUM=length(IND_NC);
C_CURR_M = PSTH_RETURN_n(Spikes.S(IND_NC,:),Infos(IND_NC,11),Start_M,End_M,SIG);



for k=1:size(W_CURR_M,2)
    Current_p_M(k,1) = ttest_NN(W_CURR_M(:,k),C_CURR_M(:,k));
end



IND_1 = find(Current_p_M<0.05);
time = Start_M:End_M;
IND_time = time(IND_1);
IND_time=IND_time(find(IND_time>175));

% find the first point that is <0.05 and is below 0.05 for at least the next 25 ms
flag=0;
for n=1:length(IND_time)
    if flag==0
        if IND_time(n+25)==IND_time(n)+25;
            Latency_ttest = IND_time(n)
            N = n;
            flag=1;
        end
    end
end



WRG_COLOUR = [255 0 0]/255;
COR_COLOUR = [0 0 255]/255;




%%%%%% start plotting here %%%%%%%%%%%%

F = figure();
subplot(2,2,1)

% Curr Wrng
IND = find(ISWRNG_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
W_CURR_NUM=length(IND);
NW_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target

% Curr Corr
IND = find(ISCORR_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
C_CURR_NUM=length(IND);
NC_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target

set(gca,'fontsize',7); ylabel('sp/s'); xlabel('time (ms)'); xlim([StartLIM_M EndLIM_M]);

hold on;
plot([Latency_ttest,Latency_ttest],ylim,'-k')
YLIM=ylim;
text(Latency_ttest,YLIM(1)+5,num2str(Latency_ttest));



subplot(2,2,3)
plot(Start_M:End_M,Current_p_M);
xlim([StartLIM_M EndLIM_M]);
box off;
hold on;
plot(xlim,[0.05 0.05],'-k')
plot([Latency_ttest Latency_ttest],ylim,'-k')
ylabel('p value from ttest');



%%%%%% to find out cell's preference
time = Start_M:End_M;
IND_new = find(time==Latency_ttest);

if nanmean(NW_M(1,IND_new+50:IND_new+100))>nanmean(NC_M(1,IND_new+50:IND_new+100))
    Pref = 'W'
else
    Pref='C'
end






%%%%%%% THEN DO ROC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



W_CURR_M; C_CURR_M;

for kk=1:length(W_CURR_M)
%     kk
    Vector1 = W_CURR_M(:,kk);
    Vector2 = C_CURR_M(:,kk);
    AUC(kk)= ROC_n(Vector1,Vector2);
    clear Vector1 Vector2;
end




IND_2 = find(AUC>0.75);
time = Start_M:End_M;
IND_time = time(IND_2);
IND_time=IND_time(find(IND_time>175));

% find the first point that is <0.05 and is below 0.05 for at least the next 25 ms
flag=0;
for n=1:length(IND_time)
    if flag==0
        if IND_time(n+25)==IND_time(n)+25;
            Latency_roc = IND_time(n)
            N = n;
            flag=1;
        end
    end
end






subplot(2,2,2)

% Curr Wrng
IND = find(ISWRNG_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
W_CURR_NUM=length(IND);
NW_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,WRG_COLOUR,1,0); % Aligned to target

% Curr Corr
IND = find(ISCORR_CURR(CHANGE:CHANGE+NUM,1)==1)+CHANGE-1;
C_CURR_NUM=length(IND);
NC_M = PSTHe_n(Spikes.S(IND,:),Infos(IND,11),Start_M,End_M,SIG,COR_COLOUR,1,0); % Aligned to target

set(gca,'fontsize',7); ylabel('sp/s'); xlabel('time (ms)'); xlim([StartLIM_M EndLIM_M]);
hold on;
plot([Latency_roc,Latency_roc],ylim,'-k')
YLIM=ylim;
text(Latency_roc,YLIM(1)+5,num2str(Latency_roc));


subplot(2,2,4)
plot(Start_M:End_M,AUC);
xlim([StartLIM_M EndLIM_M]);
box off;
hold on;
ylabel('area under ROC curve');
plot(xlim,[0.75 0.75],'-k');

plot([Latency_roc,Latency_roc],ylim,'-k')







filename = strcat(NOME,'Reward_Latency');


cd(Results_dir)
print(F, '-dpdf', filename, '-r400')

cd('C:\NAVEEN_Work\Cerebellum\Workspace\TWO');
print(F, '-dpdf', filename, '-r400')






save(POP_file,'Latency_ttest','Latency_roc','Pref','AUC','-append');
save(MERGE_file,'Latency_ttest','Latency_roc','Pref','AUC','-append');
save(ALLCELLS_file,'Latency_ttest','Latency_roc','Pref','AUC','-append');





% end