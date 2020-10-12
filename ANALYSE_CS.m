%% written by naveen at JLG on 03/17/19

%% ANALYSE_CS



%% GENERAL HISTOGRAM


BIN = 50;

F = figure();

TRIALS = 1:CHANGE;

subplot(3,3,1)
HISTOGRAM_n(Spikes.C(TRIALS),Infos(TRIALS,4),-400,400,BIN,[0 0 0],0.4,1);
hold on;
ylim([0 10])

subplot(3,3,2)
HISTOGRAM_n(Spikes.C(TRIALS),Infos(TRIALS,11),-200,600,BIN,[0 0 0],0.4,1);
hold on;
ylim([0 10])

TRIALSc = find(Infos(CHANGE:CHANGE+50,10)==1)+CHANGE-1;
TRIALSw = find(Infos(CHANGE:CHANGE+50,10)==2)+CHANGE-1;

subplot(3,3,4)
HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,4),-400,400,BIN,[0 0 1],0.4,1);
HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,4),-400,400,BIN,[1 0 0],0.4,1);
hold on;
ylim([0 10])

subplot(3,3,5)
HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,11),-200,600,BIN,[0 0 1],0.4,1);
HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,11),-200,600,BIN,[1 0 0],0.4,1);
hold on;
ylim([0 10])

temp = Infos;
Infos(1:end-1,10) = Infos(2:end,10);
% Infos(2:end,10) = Infos(1:end-1,10);

TRIALSc = find(Infos(CHANGE:CHANGE+50,10)==1)+CHANGE-1;
TRIALSw = find(Infos(CHANGE:CHANGE+50,10)==2)+CHANGE-1;

subplot(3,3,7)
HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,4),-400,400,BIN,[0 0 1],0.4,1);
HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,4),-400,400,BIN,[1 0 0],0.4,1);
hold on;
ylim([0 10])

subplot(3,3,8)
HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,11),-200,600,BIN,[0 0 1],0.4,1);
HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,11),-200,600,BIN,[1 0 0],0.4,1);
hold on;
ylim([0 10])

Infos=temp;



cd(Results_dir)
filename = 'CS_HIST';
print(F, '-dpdf', filename, '-r400')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WHAT HAPPENS AFTER REWARD

TRIALSc = find(Infos(CHANGE:CHANGE+50,10)==1)+CHANGE-1;
TRIALSw = find(Infos(CHANGE:CHANGE+50,10)==2)+CHANGE-1;

[~,HIST_C] = HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,11),100,500,BIN,[0 0 1],0.4,1);
[~,HIST_W] = HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,11),100,500,BIN,[1 0 0],0.4,1);

HIST_C_REWARD = nanmean(HIST_C);
HIST_W_REWARD = nanmean(HIST_W);


TRIALS = find(Infos(1:CHANGE,10)==1)+CHANGE-1;
[~,HIST_OT] = HISTOGRAM_n(Spikes.C(TRIALS),Infos(TRIALS,11),100,500,BIN,[1 0 0],0.4,1);
HIST_OT_REWARD = nanmean(HIST_OT);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WHAT HAPPENS AROUND DELTA

temp = find(DELTA.EpochFlag==1);
HIST_C_DELTA = NaN(3,4);
HIST_W_DELTA = NaN(3,4);
HIST_OT_DELTA = NaN(3,4);

for i=1:length(temp)
    
    TEMP = Infos;
    
    if temp(i)<=2
        align = 4;
        START = -600; END = 500;
        Infos(1:end-1,10) = Infos(2:end,10);
    end
    
    if temp(i)==3
        align = 11;
        START = -500; END = 900;
        Infos(1:end-1,10) = Infos(2:end,10);
    end
    
    if temp(i)==4
        align = 11;
        START = -500; END = 900;
    end
    
%     TTIME = START:END;
%     INDD = find(TTIME>=DELTA.START(temp(i)) & TTIME<=DELTA.END(temp(i)));
    
    TRIALSc = find(Infos(CHANGE:CHANGE+30,10)==1)+CHANGE-1;
    TRIALSw = find(Infos(CHANGE:CHANGE+30,10)==2)+CHANGE-1;
    
    [XX,HIST_C] = HISTOGRAM_n(Spikes.C(TRIALSc),Infos(TRIALSc,align),START,END,BIN,[0 0 1],0.4,0);
    [XX,HIST_W] = HISTOGRAM_n(Spikes.C(TRIALSw),Infos(TRIALSw,align),START,END,BIN,[1 0 0],0.4,0);
    
    TRIALS = find(Infos(1:CHANGE,10)==2)+CHANGE-1;
    [XX,HIST_OT] = HISTOGRAM_n(Spikes.C(TRIALS),Infos(TRIALS,align),START,END,BIN,[1 0 0],0.4,0);
    
    
    %%%% BEFORE
    INDD = find(XX>=DELTA.START(temp(i))-150 & XX<=DELTA.END(temp(i))-150);
    HIST_C_DELTA(1,temp(i)) = nanmean(HIST_C(INDD));
    HIST_W_DELTA(1,temp(i)) = nanmean(HIST_W(INDD));
    HIST_OT_DELTA(1,temp(i)) = nanmean(HIST_OT(INDD));
    
    %%%% DURING
    INDD = find(XX>=DELTA.START(temp(i)) & XX<=DELTA.END(temp(i)));
    HIST_C_DELTA(2,temp(i)) = nanmean(HIST_C(INDD));
    HIST_W_DELTA(2,temp(i)) = nanmean(HIST_W(INDD));
    HIST_OT_DELTA(2,temp(i)) = nanmean(HIST_OT(INDD));
    
    %%%% AFTER
    INDD = find(XX>=DELTA.START(temp(i))+150 & XX<=DELTA.END(temp(i))+150);
    HIST_C_DELTA(3,temp(i)) = nanmean(HIST_C(INDD));
    HIST_W_DELTA(3,temp(i)) = nanmean(HIST_W(INDD));
    HIST_OT_DELTA(3,temp(i)) = nanmean(HIST_OT(INDD));
    
    Infos=TEMP;
end


CSstats.HIST_C_DELTA = HIST_C_DELTA;
CSstats.HIST_W_DELTA = HIST_W_DELTA;
CSstats.HIST_OT_DELTA = HIST_OT_DELTA;

CSstats.HIST_C_REWARD = HIST_C_REWARD;
CSstats.HIST_W_REWARD = HIST_W_REWARD;
CSstats.HIST_OT_REWARD = HIST_OT_REWARD;


save(MERGE_file,'CSstats','-append');
save(ALLCELLS_file,'CSstats','-append');
save(POP_file,'CSstats','-append');

