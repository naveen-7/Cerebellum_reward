% written by naveen at cumc on 10/23/17

%% function ZERO_Washout



Start_T = -450;
End_T  =1050;
StartLIM_T = -400; EndLIM_T = 800;

Start_M = -750;
End_M  =750;
StartLIM_M = -500; EndLIM_M = 700;



num_trials = 15;
RANGE = CHANGE-num_trials:CHANGE+num_trials;
Start_time = Start_M;
End_time = End_M;
Sigma = 100;



% REAL STUFF ----------------------------------------------------------

Signal = Spikes.S(RANGE);
Align_time = Infos(RANGE,11);
PSTH_real = PSTH_RETURN_n(Signal,Align_time,Start_time,End_time,Sigma);
    
for i=1:500
    bef_ind = randi(num_trials,[1 5]);
    aft_ind = randi([num_trials+1 num_trials*2],[1 5]);
    
    PSTH_real_bef(i,:) = nanmean(PSTH_real(bef_ind,:));
    PSTH_real_aft(i,:) = nanmean(PSTH_real(aft_ind,:));
end

PSTH_real_bef_all = nanmean(PSTH_real_bef,1);
PSTH_real_aft_all = nanmean(PSTH_real_aft,1);


REAL = dtw(PSTH_real_bef_all,PSTH_real_aft_all);

BEF_real = trapz(PSTH_real_bef_all);
AFT_real = trapz(PSTH_real_aft_all);
REAL_area = abs(BEF_real-AFT_real);

REAL_diff = pdist2(PSTH_real_bef_all,PSTH_real_aft_all,'euclidean');

[~,KS_p] = kstest2(PSTH_real_bef_all,PSTH_real_aft_all)



% SHUFFLE STUFF -------------------------------------------------------

for i=1:5
    Signal = Spikes.S(RANGE);
    Align_time = Infos(RANGE,11);
    p = randperm(length(Signal));
    Signal = Signal(p);
    Align_time = Align_time(p);
    
    PSTH_shuffle = PSTH_RETURN_n(Signal,Align_time,Start_time,End_time,Sigma);
    
    bef_ind = randi(num_trials,[1 5]);
    aft_ind = randi([num_trials+1 num_trials*2],[1 5]);
      
    PSTH_shuffle_bef(i,:) = nanmean(PSTH_shuffle(bef_ind,:));
    PSTH_shuffle_aft(i,:) = nanmean(PSTH_shuffle(aft_ind,:));
end

PSTH_shuffle_bef_all = nanmean(PSTH_shuffle_bef,1);
PSTH_shuffle_aft_all = nanmean(PSTH_shuffle_aft,1);
clear PSTH_shuffle_bef PSTH_shuffle_aft

SHUFFLE = dtw(PSTH_shuffle_bef_all,PSTH_shuffle_aft_all);

BEF = trapz(PSTH_shuffle_bef_all);
AFT = trapz(PSTH_shuffle_aft_all);
SHUFFLE_area = abs(BEF-AFT);

SHUFFLE_diff = pdist2(PSTH_shuffle_bef_all,PSTH_shuffle_aft_all,'euclidean');



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
filename = strcat(NOME,'_Association_Learning');
print(F, '-dpdf', filename, '-r400')



%%%%%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ass_Learn.ORIGNAL = ORIGINAL;
Ass_Learn.REAL = REAL;
Ass_Learn.SHUFFLE = SHUFFLE;


Ass_Learn.REALarea = REAL_area;
Ass_Learn.SHUFFLEarea = SHUFFLE_area;

Ass_Learn.BEF_real = BEF_real;
Ass_Learn.AFT_real = AFT_real;

Ass_Learn.REAL_diff = REAL_diff;
Ass_Learn.SHUFFLE_diff = SHUFFLE_diff;


Ass_Learn.KS_p = KS_p;


save(POP_file,'Ass_Learn','-append');
save(MERGE_file,'Ass_Learn','-append');
save(ALLCELLS_file,'Ass_Learn','-append');














%% end