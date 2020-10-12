%% Function to analyse hand movement for a single session for all trials
% written by naveen at cumc on 8/29/17

function ANALYSE_HAND_MOVEMENT_POP



clc;
clear;
close all;


codes_dir = 'C:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
data_dir  = 'C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES';


disp('!!!!!  ANALYSE_HAND_MOVEMENT has started running  !!!!!')
cd(data_dir)

% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE XML DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATA_file = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATA_file)
cd(PathName1);



NOME = FileName1(1:12);


cd('C:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS')
mkdir(strcat(NOME));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS\',NOME);



SIZE = size(HAND.F,2);




%%%%%% GROSS TOTOAL HAND MOVEMENT ---------------------



F = figure();
map = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
hold on;

for i=1:SIZE
    
    clear track_temp;
    track_temp(:,2) = HAND.H{1,i};
    track_temp(:,3) = HAND.V{1,i};
    
    X = track_temp(:,3); Y = track_temp(:,2);
    colour_line_n(X,Y,0.1,parula)
    
end

title(NOME)

axis off;
colorbar
filename = strcat('ALL_',NOME);
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');

clear F;



%% barney vertical inversion
if strcmp(upper(NOME(1)),'B')
    for i=1:SIZE
        HAND.V{:,i} = nanmax(HAND.V{:,i})-HAND.V{:,i};
    end
end
%% silas vertical shift
if strcmp(upper(NOME(1)),'S')
    for i=1:SIZE
        HAND.V{:,i} = -nanmin(HAND.V{:,i})+HAND.V{:,i};
    end
end
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:SIZE
    clear TEMP1 TEMP2  track_temp
    
    track_temp(:,2) = smooth(smooth(smooth(smooth(HAND.H{1,i}))));
    track_temp(:,3) = smooth(smooth(smooth(smooth(HAND.V{1,i}))));
    
    TEMP1 = track_temp(:,2);
    TEMP1_temp = TEMP1;
    D = abs(diff(TEMP1_temp));
    [~,MAXPOINT] = nanmax(TEMP1_temp);
    D_befmax = D(1:MAXPOINT-30);
    LEVEL =nanmean(diff(TEMP1_temp))+4*nanstd(diff(TEMP1_temp));
    [L, numberOfPeaks] = bwlabel(D_befmax>LEVEL);
    if numberOfPeaks~=0
        TEMP1_temp = TEMP1_temp(find(L==numberOfPeaks)+1:end);
    end
    clear D;
    D = abs(diff(TEMP1_temp));
    temp_mean = nanmean(D(find(L==numberOfPeaks)+15:end))
    temp_std = nanstd(D(find(L==numberOfPeaks)+15:end))
    temp_temp = find(D<temp_mean+temp_std,1);
    if ~isempty(temp_temp)
        TEMP1_temp = TEMP1_temp(temp_temp:end);
    end
    temp_TEMP1_temp = max(TEMP1_temp)-TEMP1_temp;
    [~,XPOS] = findpeaks(temp_TEMP1_temp);
    [~,MINPOINT] = nanmin(TEMP1_temp);
    [~,MAXPOINT] = nanmax(TEMP1_temp);
    if ~isempty(XPOS)
        MINPOINT = XPOS(1);
    end
    if MINPOINT< MAXPOINT
        TEMP1_temp = TEMP1_temp(MINPOINT:end);
    end
    
    
    
    TEMP2 = track_temp(:,3);
    TEMP2_temp = TEMP2;
    clear D;
    D = abs(diff(TEMP2_temp));
    [~,MAXPOINT] = nanmax(TEMP2_temp);
    D_befmax = D(1:MAXPOINT-30);
    LEVEL =nanmean(diff(TEMP2_temp))+4*nanstd(diff(TEMP2_temp));
    clear L;
    [L, numberOfPeaks] = bwlabel(D_befmax>LEVEL);
    if numberOfPeaks~=0
        TEMP2_temp = TEMP2_temp(find(L==numberOfPeaks)+1:end);
    end
    if ~strcmp(upper(NOME(1)),'S')
        clear D;
        D = abs(diff(TEMP2_temp));
        temp_mean = nanmean(D(find(L==numberOfPeaks)+25:end));
        temp_std = nanstd(D(find(L==numberOfPeaks)+25:end));
        temp_temp = find(D<temp_mean+temp_std,1);
        if ~isempty(temp_temp)
            TEMP2_temp = TEMP2_temp(temp_temp:end);
        end
    end
    temp_TEMP2_temp = max(TEMP2_temp)-TEMP2_temp;
    [~,XPOS] = findpeaks(temp_TEMP2_temp);
    [~,MINPOINT] = nanmin(TEMP2_temp);
    [~,MAXPOINT] = nanmax(TEMP2_temp);
    LEVEL2 = nanmean(temp_TEMP2_temp(1:10));
     [L, numberOfPeaks] = bwlabel(temp_TEMP2_temp(1:MAXPOINT)>LEVEL2);
     
    
    if (numberOfPeaks~=0)
        MINPOINT = find(L==numberOfPeaks);
    end
    if MINPOINT< MAXPOINT
        TEMP2_temp = TEMP2_temp(MINPOINT:end);
    end
    
    
    
    
    LEN1 = length(TEMP1_temp);
    LEN2 = length(TEMP2_temp);
    MINLEN = nanmin(LEN1,LEN2);
    [MAXLEN,POSLEN] = nanmax([LEN1;LEN2]);
    diflen = MAXLEN-MINLEN;
    if POSLEN==1
        TEMP1_temp = TEMP1_temp(diflen+1:end);
    end
    if POSLEN==2
        TEMP2_temp = TEMP2_temp(diflen+1:end);
    end
    
    
    IND = diff(TEMP1_temp(1:30))==0;
    TEMP1_temp(IND) = [];
    TEMP2_temp(IND) = [];
    
    IND = diff(TEMP2_temp(1:30))==0;
    TEMP1_temp(IND) = [];
    TEMP2_temp(IND) = [];
    
%     
%     if i==5 
%         TEMP1_temp = TEMP1_temp(38:end);
%         TEMP2_temp = TEMP2_temp(38:end);
%     end
%     
%      if i==6 
%         TEMP1_temp = TEMP1(38:end);
%         TEMP2_temp = TEMP2(38:end);
%     end
    
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();
    subplot(2,2,1)
    hold on;
    plot(TEMP1); plot(TEMP1_temp);
    
    subplot(2,2,2)
    hold on;
    plot(TEMP2); plot(TEMP2_temp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    TEMP1_temp = TEMP1_temp-TEMP1_temp(1);
    TEMP2_temp = TEMP2_temp-TEMP2_temp(1);
    
    HAND_new.H{1,i} = TEMP1_temp;
    HAND_new.V{1,i} = TEMP2_temp;
    
    
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












for i=1:SIZE
    TRACK_LENGTH(i,1) = length(HAND_new.H{1,i});
end
MAX_LEN = nanmax(TRACK_LENGTH);
TRACKS=NaN(MAX_LEN,SIZE*2);

FF = figure();

suptitle(NOME)

for i=1:SIZE
    
    clear track_temp;
    track_temp(:,1) = 1:length(HAND_new.H{1,i});
    track_temp(:,2) = HAND_new.H{1,i};
    track_temp(:,3) = HAND_new.V{1,i};
    
    TRACKS(1:length(track_temp(:,2)),2*i-1) = track_temp(:,2);
    TRACKS(1:length(track_temp(:,3)),2*i) = track_temp(:,3);
    
    ax1 = subplot(2,2,1);
    plot(track_temp(:,1),track_temp(:,2))
    box off;
    YLIM1 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
    ax2 = subplot(2,2,3);
    plot(track_temp(:,1),track_temp(:,3))
    box off;
    YLIM2 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
  
    
    if i==SIZE
        % 200 fms
        linkaxes([ax1,ax2],'xy')
        subplot(2,2,1)
        axis off;
        subplot(2,2,3)
        axis off;
        hold on;
        XLIM = xlim;
        plot([XLIM(2)-50 XLIM(2)],[nanmin([YLIM1(1) YLIM2(1)]) nanmin([YLIM1(1) YLIM2(1)])],'-k')
        text(XLIM(2)-40,nanmin([YLIM1(1) YLIM2(1)])-20,'50 ms','fontsize',8);
    end
    
    XLIM_ALL(i,:)=xlim;
    
    %%% velocity
    
    
    ax3 = subplot(2,2,2);
    plot(track_temp(1:length(track_temp(:,1))-1,1),diff(track_temp(:,2)))
    box off;
    YLIM1 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
    ax4 = subplot(2,2,4);
    plot(track_temp(1:length(track_temp(:,1))-1,1),diff(track_temp(:,3)))
    box off;
    YLIM2 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
    
    
    if i==SIZE
        % 200 fms
        linkaxes([ax3,ax4],'xy')
        subplot(2,2,2)
        axis off;
        subplot(2,2,4)
        axis off;
        hold on;
        XLIM = xlim;
        plot([XLIM(2)-50 XLIM(2)],[nanmin([YLIM1(1) YLIM2(1)]) nanmin([YLIM1(1) YLIM2(1)])],'-k')
        text(XLIM(2)-40,nanmin([YLIM1(1) YLIM2(1)])-20,'50 ms','fontsize',8);
    end
    
end


for i=1:4
    subplot(2,2,i)
   xlim([ nanmin(XLIM_ALL(:,1)) nanmax(XLIM_ALL(:,2))])
end





filename = strcat('ALL_',NOME,'_HV');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400');








FF = figure();

suptitle(NOME)

H1 = subplot(2,2,1);
H_MEAN_n = nanmean(TRACKS(:,1:2:end)');
H_STD = nanstd(TRACKS(:,1:2:end)')/sqrt(SIZE);
errorline_n(1:length(H_MEAN_n),H_MEAN_n,H_STD)

H2 = subplot(2,2,3);
V_MEAN_n = nanmean(TRACKS(:,2:2:end)');
V_STD = nanstd(TRACKS(:,2:2:end)')/sqrt(SIZE);
errorline_n(1:length(V_MEAN_n),V_MEAN_n,H_STD)

linkaxes([H1, H2],'xy');

filename = strcat('MEAN_',NOME,'_HV');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400');









FILE = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES\',NOME);


try
    save(FILE,'HAND_new','H_MEAN_n','V_MEAN_n','-append');
catch
    save(FILE,'HAND_new','H_MEAN_n','V_MEAN_n');
end








% %%%%%% GROSS TOTOAL HAND MOVEMENT ---------------------
% 
% 
% 
% F = figure();
% map = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
% hold on;
% 
% for i=1:SIZE
%     
%     clear track_temp;
%     track_temp(:,2) = HAND_new.H{1,i};
%     track_temp(:,3) = HAND_new.V{1,i};
%     
%     X = track_temp(:,3); Y = track_temp(:,2);
%     colour_line_n(X,Y,0.1,parula)
%     
% end
% 
% title(NOME)
% 
% axis off;
% colorbar
% filename = strcat('ALL_',NOME);
% cd(Results_dir);
% print(F, '-dpdf', filename, '-r400');
% 
% clear F;




end
