%% LEFT-RIGHT

TRIALS = 1:CHANGE;
IND_L = find(Infos(TRIALS,9)==1);
IND_R = find(Infos(TRIALS,9)==2);

Sigma = 20;

F = figure();
subplot(2,2,1); hold on;
Signal = Spikes.S(IND_L); Align_time = Infos(IND_L,4); Start_time = -900; End_time = 200;  Colour = [1 0 0];
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,Colour,1,0);
Signal = Spikes.S(IND_R); Align_time = Infos(IND_R,4); Start_time = -900; End_time = 200;  Colour = [0 0 1];
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,Colour,1,0);

subplot(2,2,2); hold on;
Signal = Spikes.S(IND_L); Align_time = Infos(IND_L,11); Start_time = -600; End_time = 700; Colour = [1 0 0];
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,Colour,1,0);
Signal = Spikes.S(IND_R); Align_time = Infos(IND_R,11); Start_time = -600; End_time = 700;  Colour = [0 0 1];
PSTHe_n(Signal,Align_time,Start_time,End_time,Sigma,Colour,1,0);


suptitle(strcat('PSTH-SS-',NOME(1:8),'-',NOME(10:12),'-LeftRight'));

cd(Results_dir)
filename = 'PSTH_Simple spike_LR';
print(F, '-dpdf', filename, '-r400')
