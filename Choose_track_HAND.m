% Function to choose the correct track form a mixture of track outputs
% for silas hand data
% written by naveen at cumc on 11/7/17


% function Choose_track_HAND

clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('C:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW');
data_dir  = fullfile('C:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS');

cd(data_dir)



%% ONE: LOADING ALL THE FILES --------------------------------------------

% LOADING xcell FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD EXCEL FILE')
[FileName1,PathName1] = uigetfile('*.csv','File to Append');   % Open standard dialog box for retrieving files
File = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));

NOME = FileName1(1:12);
TRIAL = FileName1(14);


cd('C:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS')
mkdir(strcat(NOME));
Results_dir = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS\',NOME);




cd(PathName1);

[TEMP] = xlsread(File);

TEMP = TEMP(:,2:size(TEMP,2));
NUM_TRACKS = (size(TEMP,2)+1)/3;


TEMP = reshape(TEMP(~isnan(TEMP)),[],2*NUM_TRACKS);
clear Rotation_mat;
Rotation_mat = ([0 -1; 1 0]);
% Rotation_mat = ([1 0; 0 1]);

for i=1:NUM_TRACKS
    XY=TEMP(:,2*i-1:2*i);
    TEMP_rotated(:,2*i-1:2*i) = TEMP(:,2*i-1:2*i)*Rotation_mat;
end


figure

for i=1:NUM_TRACKS
    XY=TEMP_rotated(:,2*i-1:2*i);
    
    h1= subplot(2,2,1);
    hold on;
    plot(smooth(smooth(XY(:,1))));
    text(length(XY(:,1)),XY(end,1),num2str(i))
    
    h2=subplot(2,2,3);
    hold on;
    plot(smooth(smooth(XY(:,2))));
    text(length(XY(:,2)),XY(end,2),num2str(i))
    
end

linkaxes([h1,h2],'xy')



disp('Enter the track')
TRACK = upper(input('Here: '));

track_temp(:,1) = 1:length(TEMP_rotated);
track_temp(:,2:3) = TEMP_rotated(:,2*TRACK-1:2*TRACK);
track_temp(:,2)= smooth(smooth(track_temp(:,2)));
track_temp(:,3)= smooth(smooth(track_temp(:,3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



F = figure();
% plot((track_temp(:,3)),(track_temp(:,2)))
map = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
X = (track_temp(:,2)); Y = (track_temp(:,3));
colour_line_n(X,Y,3,parula)
axis off;
colorbar
filename = strcat(NOME,'_',num2str(TRIAL));
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');
clear F;





F = figure();
X = track_temp(:,2); Y = track_temp(:,3);
plot(((track_temp(:,2))),((track_temp(:,3))),'-r','linewidth',2)
axis off;
filename = strcat(NOME,'_',num2str(TRIAL),'_SAMECOLOR');
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');
clear F;





FF = figure();

subplot(2,2,1)
plot(track_temp(:,1),track_temp(:,2))
box off;
YLIM1 = ylim;
xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
axis off;

subplot(2,2,3)
plot(track_temp(:,1),track_temp(:,3))
box off;
YLIM2 = ylim;
xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
axis off;


subplot(2,2,1)
ylim([nanmin([YLIM1(1) YLIM2(1)]) nanmax([YLIM1(2) YLIM2(2)]) ]);
subplot(2,2,3)
ylim([nanmin([YLIM1(1) YLIM2(1)]) nanmax([YLIM1(2) YLIM2(2)]) ]);

% 250 fms
hold on;
XLIM = xlim;
plot([XLIM(2)-50 XLIM(2)],[nanmin([YLIM1(1) YLIM2(1)]) nanmin([YLIM1(1) YLIM2(1)])],'-k')
text(XLIM(2)-40,nanmin([YLIM1(1) YLIM2(1)])-20,'50 ms','fontsize',8);




%%% velocity


subplot(2,2,2)
plot(track_temp(1:length(track_temp(:,1))-1,1),diff(track_temp(:,2)),'color',[217 83 25]/255)
box off;
YLIM1 = ylim;
xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
axis off;

subplot(2,2,4)
plot(track_temp(1:length(track_temp(:,1))-1,1),diff(track_temp(:,3)),'color',[217 83 25]/255)
box off;
YLIM2 = ylim;
xlim([track_temp(1,1) track_temp(length(track_temp),1)]);
axis off;



subplot(2,2,2)
ylim([nanmin([YLIM1(1) YLIM2(1)]) nanmax([YLIM1(2) YLIM2(2)]) ]);
subplot(2,2,4)
ylim([nanmin([YLIM1(1) YLIM2(1)]) nanmax([YLIM1(2) YLIM2(2)]) ]);


% 200 fms
hold on;
XLIM = xlim;
plot([XLIM(2)-50 XLIM(2)],[nanmin([YLIM1(1) YLIM2(1)]) nanmin([YLIM1(1) YLIM2(1)])],'-k')
text(XLIM(2)-40,nanmin([YLIM1(1) YLIM2(1)])-20,'50 ms','fontsize',8);



filename = strcat(NOME,'_',num2str(TRIAL),'_HV');
cd(Results_dir);
print(FF, '-dpdf', filename, '-r400');


FILE = strcat('C:\NAVEEN_Work\Cerebellum\Data\videos\ALL_FILES\',NOME);

try
load(FILE)
catch
    1;
end


HAND.F{str2num(TRIAL)} = track_temp(:,1);
HAND.H{str2num(TRIAL)} = track_temp(:,2);
HAND.V{str2num(TRIAL)} = track_temp(:,3);

VEL.F{str2num(TRIAL)} = track_temp(1:length(track_temp(:,1))-1,1);
VEL.H{str2num(TRIAL)} = diff(track_temp(:,2));
VEL.V{str2num(TRIAL)} = diff(track_temp(:,3));


try
    save(FILE,'HAND','VEL','-append');
catch
    save(FILE,'HAND','VEL');
end





% end

