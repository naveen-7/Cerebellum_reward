
%% function ZERO_motor_learning

% Written by naveen on 10/11/17 at cumc

% Take 10 trials before change and 10 trials after change
% take 5 random trials in bef and 5 random trials in aft and take mean
% compute DTW for the mean
% repeat 500 times
% Randomize the trials and repeat





% Create a structure
% bars then dowels

clear spikes;
spikes.S(1:length(Spikes.B),1) = Spikes.B;
spikes.S(length(Spikes.B)+1:length(Spikes.B)+length(Spikes.D),1) = Spikes.D;

CHANGE=11;

clear infos
infos(1:size(INFOS.B,1),:) = INFOS.B;
infos(size(INFOS.B,1)+1:size(INFOS.B,1)+size(INFOS.D,1),:) = INFOS.D;










Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;



num_trials = 10;
RANGE = CHANGE-num_trials:CHANGE+num_trials-1;
Start_time = Start_M;
End_time = End_M;
Sigma = 30;



% REAL STUFF ----------------------------------------------------------

Signal = spikes.S(RANGE);
Align_time = infos(RANGE,11);
PSTH_real = PSTH_RETURN_n(Signal,Align_time,Start_time,End_time,Sigma);
    
for i=1:500
    bef_ind = randi(num_trials,[1 5]);
    aft_ind = randi([num_trials+1 num_trials*2],[1 5]);
    
    PSTH_real_bef(i,:) = nanmean(PSTH_real(bef_ind,:));
    PSTH_real_aft(i,:) = nanmean(PSTH_real(aft_ind,:));
end

PSTH_real_bef_all = nanmean(PSTH_real_bef);
PSTH_real_aft_all = nanmean(PSTH_real_aft);

REAL = dtw(PSTH_real_bef_all,PSTH_real_aft_all);




% SHUFFLE STUFF -------------------------------------------------------

for i=1:1000
    Signal = spikes.S(RANGE);
    p = randperm(length(Signal));
    Signal = Signal(p);
    Align_time = Align_time(p);
    
    PSTH_shuffle = PSTH_RETURN_n(Signal,Align_time,Start_time,End_time,Sigma);
    
    bef_ind = randi(num_trials,[1 5]);
    aft_ind = randi([num_trials+1 num_trials*2],[1 5]);
      
    PSTH_shuffle_bef(i,:) = nanmean(PSTH_shuffle(bef_ind,:));
    PSTH_shuffle_aft(i,:) = nanmean(PSTH_shuffle(aft_ind,:));
end

PSTH_shuffle_bef_all = nanmean(PSTH_shuffle_bef);
PSTH_shuffle_aft_all = nanmean(PSTH_shuffle_aft);
clear PSTH_shuffle_bef PSTH_shuffle_aft

% SHUFFLE = dtw(PSTH_shuffle_bef_all,PSTH_shuffle_aft_all);
SHUFFLE = dtw(PSTH_shuffle_bef_all,PSTH_shuffle_aft_all);



F = figure();

time = Start_time:End_time;

subplot(2,2,1)
hold on;
PSTH_real_bef_real = nanmean(PSTH_real(1:num_trials,:));
PSTH_real_aft_real = nanmean(PSTH_real(num_trials+1:end,:));
ORIGINAL = dtw(PSTH_real_bef_real,PSTH_real_aft_real);
plot(time,PSTH_real_bef_real);
plot(time,PSTH_real_aft_real);
xlim([Start_time End_time])
YLIM = ylim;
title('original')
text(Start_time+50,YLIM(2)-10,num2str(round(ORIGINAL)))

% % % % to illustrate only
% % % for i=1:length(PSTH_real_aft_all)
% % %     PSTH_real_aft_all(i) = PSTH_real_bef_all(i)+rand;
% % % end
% % % PSTH_real_aft_all = smooth(PSTH_real_aft_all);

subplot(2,2,3)
hold on;
plot(time,PSTH_real_bef_all);
plot(time,PSTH_real_aft_all);
xlim([Start_time End_time])
ylim([YLIM(1) YLIM(2)])
title('real random')
text(Start_time+50,YLIM(2)-10,num2str(round(REAL)))

subplot(2,2,4)
hold on;
plot(time,PSTH_shuffle_bef_all);
plot(time,PSTH_shuffle_aft_all);
xlim([Start_time End_time])
ylim([YLIM(1) YLIM(2)])
title('shuffle random')
text(Start_time+50,YLIM(2)-10,num2str(round(SHUFFLE)))


cd(Results_dir)
filename = strcat(NOME,'_Motor_Learning');
print(F, '-dpdf', filename, '-r400')



%%%%%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Mot_Learn.ORIGNAL = ORIGINAL;
Mot_Learn.REAL = REAL;
Mot_Learn.SHUFFLE = SHUFFLE;

save(POP_file,'Mot_Learn','-append');
save(MERGE_file,'Mot_Learn','-append');
save(ALLCELLS_file,'Mot_Learn','-append');




% end