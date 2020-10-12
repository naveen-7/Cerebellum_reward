%% This function does something similar to MRD using a method like bootstrap
% This is to show that the delta epoch is a real thing
% Created by Naveen on 10/17/18 at JLG

function [NULL_mean,SHUFF_mean,p] = shoestrap_n(Infos,Spikes,CHANGE,DELTA,cycles,NUM)

% Infos*      : n x x : Infos structure
% Spikes*     : n x x : Spikes structure
% CHANGE*     : 1 x 1 : CHANGE from one condition to another, usualy from OT to N
% DELTA*      : x x x : DELTA structure
% cycles      : 1 x 1 : number of cycles to run [Default: 250]
% NUM         : 1 x 1 : number of trials after change to consider [Default: 20]

if nargin<4
    error('Incomplete input to the function shoestrap_n');
elseif nargin==4
    varargin{1} = Infos;
    varargin{2} = Spikes;
    varargin{3} = CHANGE;
    varargin{4} = DELTA;
    cycles      = 250;
    NUM         = 20;
elseif nargin==5
    varargin{1} = Infos;
    varargin{2} = Spikes;
    varargin{3} = CHANGE;
    varargin{4} = DELTA;
    varargin{5} = cycles;
    NUM         = 20;
elseif nargin==6
    varargin{1} = Infos;
    varargin{2} = Spikes;
    varargin{3} = CHANGE;
    varargin{4} = DELTA;
    varargin{5} = cycles;
    varargin{6} = NUM;
else
    error('Too many inputs to the function shoestrap_n');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIG = 40;

Start_T = -450;
End_T  =1250;
StartLIM_T = -400; EndLIM_T = 1000;

Start_M = -750;
End_M  =950;
StartLIM_M = -500; EndLIM_M = 900;





loc = find(DELTA.EpochFlag==1);

if length(loc)==1 %%%%%%%%%%%%%%%%% ONLY DO SINGLE DELTA neurons
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:1 initialize everything according to delta location
    if loc<=2
        Align = 4;
        Start = Start_T; End = End_T; StartLIM = StartLIM_T; EndLIM = EndLIM_T;
        psthtemp=PSTH_RETURN_n(Spikes.S(CHANGE:CHANGE+NUM-1,1),Infos(CHANGE:CHANGE+NUM-1,11),Start,End,SIG);
        time = Start:End;
        SDF_tru = psthtemp(:,find(time==StartLIM): find(time==EndLIM));
        D_srt = DELTA.START(loc);
        D_end = DELTA.END(loc);
        time_tru = StartLIM:EndLIM;
    end
    if loc==3|loc==4
        Align = 11;
        Start = Start_M; End = End_M; StartLIM = StartLIM_M; EndLIM = EndLIM_M;
        psthtemp=PSTH_RETURN_n(Spikes.S(CHANGE:CHANGE+NUM-1,1),Infos(CHANGE:CHANGE+NUM-1,11),Start,End,SIG);
        time = Start:End;
        SDF_tru = psthtemp(:,find(time==StartLIM): find(time==EndLIM));
        D_srt = DELTA.START(loc);
        D_end = DELTA.END(loc);
        time_tru = StartLIM:EndLIM;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:2 create matrix of PSTH and trial outcomes amd time
    
    if loc<=3
        VAL_tru = Infos(CHANGE-1:CHANGE+NUM,10);
    end
    
    if loc==4
        VAL_tru = Infos(CHANGE:CHANGE+NUM-1,10);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:3 do shoestrap for null data
    
    for i=1:cycles
        SDF_use = SDF_tru;
        VAL_use = VAL_tru;
        SHUFF_IND = randperm(20);
        SDF_use = SDF_use(SHUFF_IND,:);
        VAL_use = VAL_use(SHUFF_IND,:);
        
        corr_ind = find(VAL_use==1,5);
        wrng_ind = find(VAL_use==2,5);
        EPOCH = find(time_tru==D_srt):find(time_tru==D_end);
        clear COMPARE;
        COMPARE(1,:)=nanmean(SDF_use(corr_ind,EPOCH));
        COMPARE(2,:)=nanmean(SDF_use(wrng_ind,EPOCH));
        
        
        NULL(i,1) = sqrt(sum((COMPARE(1,:)-COMPARE(2,:)).^2));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:4 do shoestrap for shuffle data
    
    for i=1:cycles
        SDF_use = SDF_tru;
        VAL_use = VAL_tru;
        SHUFF_IND = randperm(20);
        %SDF_use = SDF_use(SHUFF_IND,:);
        VAL_use = VAL_use(SHUFF_IND,:);
        
        corr_ind = find(VAL_use==1,5);
        wrng_ind = find(VAL_use==2,5);
        EPOCH = find(time_tru==D_srt):find(time_tru==D_end);
        clear COMPARE;
        COMPARE(1,:)=nanmean(SDF_use(corr_ind,EPOCH));
        COMPARE(2,:)=nanmean(SDF_use(wrng_ind,EPOCH));
        
        
        SHUFF(i,1) = sqrt(sum((COMPARE(1,:)-COMPARE(2,:)).^2));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:4 plotting nd stuff
    
    F = figure()
    subplot(2,2,1)
    hold on;
    NULL_color = [209 173 21]/255;
    SHUFF_color = [102 102 102]/255;
    BIN_W = 40;
    h = histogram(NULL,'BinWidth',BIN_W);
    h.FaceColor = NULL_color;
    h = histogram(SHUFF,'BinWidth',BIN_W);
    h.FaceColor = SHUFF_color;
    box off;
    xlabel('rms distance(a.u.)');
    ylabel('frequency');
    YLIM = ylim;
    plot([nanmean(SHUFF) nanmean(NULL)],[YLIM(2) YLIM(2)],'-k')
    text(nanmean(SHUFF)+0.25*(nanmean(NULL)-nanmean(SHUFF)), YLIM(2), star_n(p));
    
    
    cd(Results_dir)
    filename = strcat(NOME,'DELTA_MRD');
    print(F, '-dpdf', filename, '-r400')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP:5 retrun these stuff
    NULL_mean = nanmean(NULL);
    SHUFF_mean = nanmean(SHUFF);
    p = stats_test_n(NULL,SHUFF);
    
    
    
end



end
