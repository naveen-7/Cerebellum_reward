% written by naveen at cumc on 12/22/17

%% function ZERO_Washout


% Take 10 trials before change and 10 trials after change
% take 5 random trials in bef and 5 random trials in aft and take mean
% compute distance for the mean
% repeat 500 times
% calculate the 'maximum allowed variance' in the original data in the same way

NUM = 50;

Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 400;
TIME_T = Start_T:End_T;

Start_M = -750;
End_M  =750;
StartLIM_M = -200; EndLIM_M = 600;
TIME_M = Start_M:End_M;


num_trials = 20;
NN = 10;
RANGE = CHANGE-num_trials:CHANGE+num_trials-1;
Start_time = Start_M;
End_time = End_M;
Sigma = 60;



% REAL STUFF ----------------------------------------------------------

clear REAL_diff1 REAL_diff2
clear PSTH_real_bef PSTH_real_aft PSTH_real_bef_all_M PSTH_real_aft_all_M

Signal = Spikes.S(RANGE);
Align_time = Infos(RANGE,11);
PSTH_real_M = PSTH_RETURN_n(Signal,Align_time,Start_M,End_M,Sigma);
PSTH_real_M = PSTH_real_M(:,find(TIME_M==StartLIM_M):find(TIME_M==EndLIM_M));

for i=1:NUM
    TEMPPPY = randperm(num_trials); bef_ind = TEMPPPY(1:NN);
    TEMPPPY = num_trials+randperm(num_trials); aft_ind = TEMPPPY(1:NN);
    
    PSTH_real_bef(i,:) = nanmean(PSTH_real_M(bef_ind,:));
    PSTH_real_aft(i,:) = nanmean(PSTH_real_M(aft_ind,:));
    
    REAL_diff1(i,:) = pdist2(PSTH_real_bef(i,:),PSTH_real_aft(i,:),'euclidean');
end


PSTH_real_bef_all_M = nanmean(PSTH_real_bef,1);
PSTH_real_aft_all_M = nanmean(PSTH_real_aft,1);



clear PSTH_real_bef PSTH_real_aft PSTH_real_bef_all_T PSTH_real_aft_all_T
Signal = Spikes.S(RANGE);
Align_time = Infos(RANGE,4);
PSTH_real_T = PSTH_RETURN_n(Signal,Align_time,Start_T,End_T,Sigma);
PSTH_real_T = PSTH_real_T(:,find(TIME_T==StartLIM_T):find(TIME_T==EndLIM_T));
    
for i=1:NUM
    TEMPPPY = randperm(num_trials); bef_ind = TEMPPPY(1:NN);
    TEMPPPY = num_trials+randperm(num_trials); aft_ind = TEMPPPY(1:NN);
    
    PSTH_real_bef(i,:) = nanmean(PSTH_real_T(bef_ind,:));
    PSTH_real_aft(i,:) = nanmean(PSTH_real_T(aft_ind,:));
    
    REAL_diff2(i,:) = pdist2(PSTH_real_bef(i,:),PSTH_real_aft(i,:),'euclidean');
end

PSTH_real_bef_all_T = nanmean(PSTH_real_bef,1);
PSTH_real_aft_all_T = nanmean(PSTH_real_aft,1);






MIN = nanmin([nanmin(nanmin(PSTH_real_bef_all_M)) nanmin(nanmin(PSTH_real_bef_all_M))  nanmin(nanmin(PSTH_real_bef_all_T))  nanmin(nanmin(PSTH_real_bef_all_T))]);
MAX = nanmax([nanmax(nanmax(PSTH_real_bef_all_M)) nanmax(nanmax(PSTH_real_bef_all_M))  nanmax(nanmax(PSTH_real_bef_all_T))  nanmax(nanmax(PSTH_real_bef_all_T))]);


PSTH_real_bef_all_M = PSTH_real_bef_all_M-MIN;
PSTH_real_bef_all_M = PSTH_real_bef_all_M/MAX;

PSTH_real_aft_all_M = PSTH_real_aft_all_M-MIN;
PSTH_real_aft_all_M = PSTH_real_aft_all_M/MAX;

PSTH_real_bef_all_T = PSTH_real_bef_all_T-MIN;
PSTH_real_bef_all_T = PSTH_real_bef_all_T/MAX;

PSTH_real_aft_all_T = PSTH_real_aft_all_T-MIN;
PSTH_real_aft_all_T = PSTH_real_aft_all_T/MAX;



REAL_Diff1 = pdist2(PSTH_real_bef_all_M,PSTH_real_aft_all_M,'euclidean');
REAL_Diff2 = pdist2(PSTH_real_bef_all_T,PSTH_real_aft_all_T,'euclidean');




%%%%%%%%%%%%%%%%%%
REAL_diff = nansum([nanmean(REAL_diff1) nanmean(REAL_diff2)]);
REAL_Diff = nansum([nanmean(REAL_Diff1) nanmean(REAL_Diff2)]); % normalized PSTH
%%%%%%%%%%%%%%%%%%








clear PSTH_mav_bef PSTH_mav_aft MAV_diff1 MAV_diff2

% MAX ALLOWABLE VARIANCE STUFF ----------------------------------------------------------
RANGE = 1:CHANGE;
Signal = Spikes.S(RANGE);
Align_time = Infos(RANGE,11);
PSTH_mav_M = PSTH_RETURN_n(Signal,Align_time,Start_M,End_M,Sigma);
PSTH_mav_M = PSTH_mav_M(:,find(TIME_M==StartLIM_M):find(TIME_M==EndLIM_M));

for i=1:NUM

    temp_now =randperm(CHANGE);
     bef_ind = temp_now(1:NN);
     aft_ind = temp_now(NN+1:num_trials);
    
    PSTH_mav_bef(i,:) = nanmean(PSTH_mav_M(bef_ind,:));
    PSTH_mav_aft(i,:) = nanmean(PSTH_mav_M(aft_ind,:));
   
    MAV_diff1(i,:) = pdist2(PSTH_mav_bef(i,:),PSTH_mav_aft(i,:),'euclidean');
end


PSTH_mav_bef_all_M = nanmean(PSTH_mav_bef,1);
PSTH_mav_aft_all_M = nanmean(PSTH_mav_aft,1);

% MAV_Diff1 = pdist2(PSTH_mav_bef_all_M,PSTH_mav_aft_all_M,'euclidean');






clear PSTH_mav_bef PSTH_mav_aft
Signal = Spikes.S(RANGE);
Align_time = Infos(RANGE,4);
PSTH_mav_T = PSTH_RETURN_n(Signal,Align_time,Start_T,End_T,Sigma);
PSTH_mav_T = PSTH_mav_T(:,find(TIME_T==StartLIM_T):find(TIME_T==EndLIM_T));

for i=1:NUM
     temp_now =randperm(CHANGE);
     bef_ind = temp_now(1:NN);
     aft_ind = temp_now(NN+1:num_trials);
    
    PSTH_mav_bef(i,:) = nanmean(PSTH_mav_T(bef_ind,:));
    PSTH_mav_aft(i,:) = nanmean(PSTH_mav_T(aft_ind,:));
   
    MAV_diff2(i,:) = pdist2(PSTH_mav_bef(i,:),PSTH_mav_aft(i,:),'euclidean');
end


PSTH_mav_bef_all_T = nanmean(PSTH_mav_bef,1);
PSTH_mav_aft_all_T = nanmean(PSTH_mav_aft,1);

% MAV_Diff2 = pdist2(PSTH_mav_bef_all_T,PSTH_mav_aft_all_T,'euclidean');





MIN = nanmin([nanmin(nanmin(PSTH_mav_bef_all_M)) nanmin(nanmin(PSTH_mav_bef_all_M))  nanmin(nanmin(PSTH_mav_bef_all_T))  nanmin(nanmin(PSTH_mav_bef_all_T))]);
MAX = nanmax([nanmax(nanmax(PSTH_mav_bef_all_M)) nanmax(nanmax(PSTH_mav_bef_all_M))  nanmax(nanmax(PSTH_mav_bef_all_T))  nanmax(nanmax(PSTH_mav_bef_all_T))]);


PSTH_mav_bef_all_M = PSTH_mav_bef_all_M-MIN;
PSTH_mav_bef_all_M = PSTH_mav_bef_all_M/MAX;

PSTH_mav_aft_all_M = PSTH_mav_aft_all_M-MIN;
PSTH_mav_aft_all_M = PSTH_mav_aft_all_M/MAX;

PSTH_mav_bef_all_T = PSTH_mav_bef_all_T-MIN;
PSTH_mav_bef_all_T = PSTH_mav_bef_all_T/MAX;

PSTH_mav_aft_all_T = PSTH_mav_aft_all_T-MIN;
PSTH_mav_aft_all_T = PSTH_mav_aft_all_T/MAX;



MAV_Diff1 = pdist2(PSTH_mav_bef_all_M,PSTH_mav_aft_all_M,'euclidean');
MAV_Diff2 = pdist2(PSTH_mav_bef_all_T,PSTH_mav_aft_all_T,'euclidean');






%%%%%%%%%%%%%%%%%%
MAV_diff = nansum([nanmean(MAV_diff1) nanmean(MAV_diff2)]);
MAV_Diff = nansum([nanmean(MAV_Diff1) nanmean(MAV_Diff2)]);
%%%%%%%%%%%%%%%%%%







MAV_all = [MAV_diff1; MAV_diff2];
REAL_all = [REAL_diff1; REAL_diff2];

REAL_color = [33 178 128]/255;
MAV_color = [102 102 102]/255;

VAL_p = ttest_NN(MAV_all,REAL_all); % ttest p value
VAL_R = ROC_n(MAV_all,REAL_all); % ROC value


BIN_W = 20;

F = figure;
subplot(2,2,1,'box','off','tickdir','out','Linewidth',0.3,'FontSize',7)
hold on;
h = histogram(REAL_all,'BinWidth',BIN_W);
h.FaceColor = REAL_color;
h = histogram(MAV_all,'BinWidth',BIN_W);
h.FaceColor = MAV_color;
box off;
xlabel('rms distance(a.u.)');
ylabel('frequency');
YLIM = ylim;
plot([nanmean(MAV_all) nanmean(REAL_all)],[YLIM(2) YLIM(2)],'-k')
text(nanmean(MAV_all)+0.25*(nanmean(REAL_all)-nanmean(MAV_all)), YLIM(2), star_n(VAL_p));


subplot(2,2,2)
hold on;
text(0.5,0.5,strcat('p = ',num2str(VAL_p)));
text(0.5,0.4,strcat('AUC = ',num2str(VAL_R)));
axis off;
box off;

cd(Results_dir)
filename = strcat(NOME,'_Association_Learning_MAV_STATS');
print(F, '-dpdf', filename, '-r400')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


F = figure();

time = StartLIM_T:EndLIM_T;

subplot(3,3,1)
hold on;
PSTH_real_bef_real = nanmean(PSTH_real_T(1:num_trials,:));
PSTH_real_aft_real = nanmean(PSTH_real_T(num_trials+1:end,:));
ORIGINAL = pdist2(PSTH_real_bef_real,PSTH_real_aft_real,'euclidean');
plot(time,PSTH_real_bef_real);
plot(time,PSTH_real_aft_real);
xlim([time(1) time(end)])
YLIM = ylim;
title('original')
text(Start_time+700,YLIM(2)-10,num2str(round(ORIGINAL)))

subplot(3,3,4)
hold on;
plot(time,PSTH_real_bef_all_T);
plot(time,PSTH_real_aft_all_T);
xlim([time(1) time(end)])
ylim([0 1])
% ylim([YLIM(1) YLIM(2)])
title('real')
text(Start_time+700,YLIM(2)-10,num2str(round(REAL_Diff2)))

subplot(3,3,7)
hold on;
plot(time,PSTH_mav_bef_all_T);
plot(time,PSTH_mav_aft_all_T);
xlim([time(1) time(end)])
ylim([0 1])
% ylim([YLIM(1) YLIM(2)])
title('mav')
text(Start_time+700,YLIM(2)-10,num2str(round(MAV_Diff2)))




time = StartLIM_M:EndLIM_M;

subplot(3,3,2)
hold on;
PSTH_real_bef_real = nanmean(PSTH_real_M(1:num_trials,:));
PSTH_real_aft_real = nanmean(PSTH_real_M(num_trials+1:end,:));
ORIGINAL = pdist2(PSTH_real_bef_real,PSTH_real_aft_real,'euclidean');
plot(time,PSTH_real_bef_real);
plot(time,PSTH_real_aft_real);
xlim([time(1) time(end)])
YLIM = ylim;
title('original')
text(Start_time+700,YLIM(2)-10,num2str(round(ORIGINAL)))

subplot(3,3,5)
hold on;
plot(time,PSTH_real_bef_all_M);
plot(time,PSTH_real_aft_all_M);
xlim([time(1) time(end)])
ylim([0 1])
% ylim([YLIM(1) YLIM(2)])
title('real')
text(Start_time+700,YLIM(2)-10,num2str(round(REAL_Diff1)))

subplot(3,3,8)
hold on;
plot(time,PSTH_mav_bef_all_M);
plot(time,PSTH_mav_aft_all_M);
xlim([time(1) time(end)])
ylim([0 1])
% ylim([YLIM(1) YLIM(2)])
title('mav')
text(Start_time+700,YLIM(2)-10,num2str(round(MAV_Diff1)))






subplot(3,3,6)
hold on;
text(0.5,0.5,num2str(round(REAL_Diff)))
axis off;

subplot(3,3,9)
hold on;
text(0.5,0.5,num2str(round(MAV_Diff)))
axis off;



cd(Results_dir)
filename = strcat(NOME,'_Association_Learning_MAV');
print(F, '-dpdf', filename, '-r400')



%%%%%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ass_Learn.REAL_mav = REAL_Diff;
Ass_Learn.MAV_mav = MAV_Diff;


Ass_Learn.VAL_p = VAL_p;
Ass_Learn.VAL_R = VAL_R;


save(POP_file,'Ass_Learn','-append');
save(MERGE_file,'Ass_Learn','-append');
save(ALLCELLS_file,'Ass_Learn','-append');




% end