%% Function to analyse hand movement for a single session for all trials
% written by naveen at cumc on 1/4/18

function ANALYSE_HAND_MOVEMENT_POP_manual



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


FILE = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES\',NOME);


% 
% 
% %%%%% ADDITION %
% %%%%%%%%%%%%%%%%
% 
% PARENT1 = 4;
% PARENT2 = 1;
% CHILD = 5;
% 
% clear X Y
% X = cell2mat(HAND.H(1,PARENT1));
% Y = cell2mat(HAND.H(1,PARENT2));
% MINLEN = nanmin([length(X),length(Y)]);
% HAND.H(1,CHILD) = {nanmean([X(1:MINLEN),Y(1:MINLEN)],2)};
% X = cell2mat(HAND.V(1,PARENT1));
% Y = cell2mat(HAND.V(1,PARENT2));
% HAND.V(1,CHILD) = {nanmean([X(1:MINLEN),Y(1:MINLEN)],2)};
% HAND.F(1,CHILD) = {[0:MINLEN-1]'};
% 
% save(FILE,'HAND','-append');
% 
% %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%







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


if exist ('CUT_all') & exist('HAND_chosen')
    disp('"CUT_all" exists');
    disp('run anyway?');
    FLAG = input('enter 0 or 1 only: ');
    
    if FLAG==0
        disp('>>> Proceeding AUTOMATICALLY >>>');
    end
end


if ~exist ('CUT_all') | ~exist('HAND_chosen')| FLAG==1
    disp('"CUT_all" does not exist');
    disp('>>> Proceeding MANUALLY >>>');
    
    clear HAND_chosen;
    i=1;
    while i<=SIZE
        clear TEMP1 TEMP2  track_temp
        
        track_temp(:,2) = HAND.H{1,i};
        track_temp(:,3) = HAND.V{1,i};
        
        track_temp(:,2) = -nanmin(track_temp(:,2))+track_temp(:,2);
        track_temp(:,3) = -nanmin(track_temp(:,3))+track_temp(:,3);
        
        figure();
        suptitle(strcat('trial : ',num2str(i)));
        subplot(2,2,1)
        hold on;
        plot(track_temp(:,2));
        %     plot(abs(diff(track_temp(:,2))*10));
        subplot(2,2,2)
        hold on;
        plot(track_temp(:,3));
        %     plot(abs(diff(track_temp(:,3))*10));
        
        
        disp('Enter cut number');
        CUT = input('here: ');
        
        HAND_new.H{1,i} = track_temp(CUT:end,2);
        HAND_new.V{1,i} = track_temp(CUT:end,3);
        
        
        subplot(2,2,1)
        plot(HAND_new.H{1,i});
        subplot(2,2,2)
        plot(HAND_new.V{1,i});
        
        
        IN = upper(input('Proceed? ','s'));
        if strcmp(IN,'Y')
            i=i+1;
            delete(gcf);
        end
        if strcmp(IN,'N')
            1;
        end
        
        
        CUT_all(i,1)=CUT;
    end  



try
    save(FILE,'CUT_all','-append');
catch
    save(FILE,'CUT_all');
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSING 5

COLOURS = get(gca,'colororder');
COLOURS=[COLOURS;COLOURS;COLOURS];
close all;

F = figure();

for i=1:SIZE
    
    
    clear track_temp;
    track_temp(:,1) = 1:length(HAND_new.H{1,i});
    track_temp(:,2) = HAND_new.H{1,i};
    track_temp(:,3) = HAND_new.V{1,i};
    
    subplot(4,SIZE,i)
    text(0.2,0.5,num2str(i));
    axis off;
    
    ax1 = subplot(4,SIZE,i+SIZE);
    hold on;
    plot(track_temp(:,1),track_temp(:,2),'color',COLOURS(i,:));
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    box off;
    axis off;
    
    ax2 = subplot(4,SIZE,i+2*SIZE);
    hold on;
    plot(track_temp(:,1),track_temp(:,3),'color',COLOURS(i,:));
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    box off;
    axis off;
    
    linkaxes([ax1,ax2],'xy')
    
end


filename = strcat('ALL_',NOME,'_HV_seperate');
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    plot(track_temp(1:length(track_temp(:,1))-1,1),smooth(smooth(diff(track_temp(:,2)))))
    box off;
    YLIM1 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
    ax4 = subplot(2,2,4);
    plot(track_temp(1:length(track_temp(:,1))-1,1),smooth(smooth(diff(track_temp(:,3)))))
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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
disp('ENTER 5 trials seperated by a dot (.)');
T = input('here: ','s');
C = strsplit(T,'.');

for count=1:5
    num = str2num(cell2mat(C(count)))
    HAND_chosen.H{1,count} = HAND_new.H{1,num};
    HAND_chosen.V{1,count} = HAND_new.V{1,num};
end
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');




for count=1:5
    HAND_chosen.H{1,count} = HAND_chosen.H{1,count}-nanmin(HAND_chosen.H{1,count});
    HAND_chosen.V{1,count} = HAND_chosen.V{1,count}-nanmin(HAND_chosen.H{1,count});
end



for count=1:5
    HAND_new.H{1,count} = HAND_new.H{1,count}-nanmin(HAND_new.H{1,count});
    HAND_new.V{1,count} = HAND_new.V{1,count}-nanmin(HAND_new.H{1,count});
end


save(FILE,'HAND_chosen','HAND_new','-append');

end
%% %% %% %%%%%%%%%%%%%%%%%%%%


%% PRINT CHOSEN 5



SIZE = 5;

for i=1:5
    TRACK_LENGTH_ch(i,1) = length(HAND_chosen.H{1,i});
end
MAX_LEN = nanmax(TRACK_LENGTH_ch);
TRACKS=NaN(MAX_LEN,SIZE*2);


FF = figure();
suptitle(NOME)

for i=1:SIZE
    
    clear track_temp;
    track_temp(:,1) = 1:length(HAND_chosen.H{1,i});
    track_temp(:,2) = HAND_chosen.H{1,i};
    track_temp(:,3) = HAND_chosen.V{1,i};
    
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
    plot(track_temp(1:length(track_temp(:,1))-1,1),smooth(smooth(diff(track_temp(:,2)))))
    box off;
    YLIM1 = ylim;
    xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
    hold on;
    
    ax4 = subplot(2,2,4);
    plot(track_temp(1:length(track_temp(:,1))-1,1),smooth(smooth(diff(track_temp(:,3)))))
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


filename = strcat('ALL_',NOME,'_HV_CHOSEN');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400');










FF = figure();

suptitle(NOME)

H1 = subplot(2,2,1);
H_MEAN = nanmean(TRACKS(:,1:2:end)');
H_STD = nanstd(TRACKS(:,1:2:end)')/sqrt(SIZE);
errorline_n(1:length(H_MEAN),H_MEAN,H_STD)

H2 = subplot(2,2,3);
V_MEAN = nanmean(TRACKS(:,2:2:end)');
V_STD = nanstd(TRACKS(:,2:2:end)')/sqrt(SIZE);
errorline_n(1:length(V_MEAN),V_MEAN,H_STD)

linkaxes([H1, H2],'xy');

filename = strcat('MEAN_',NOME,'_HV');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400');





H_MEAN_ch = H_MEAN;
V_MEAN_ch = V_MEAN;


save(FILE,'H_MEAN_ch','V_MEAN_ch','-append');









end
