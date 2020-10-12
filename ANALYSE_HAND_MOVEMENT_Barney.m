
%% Function to analyse hand movement trial by trial
% written by naveen at cumc on 8/26/17

function ANALYSE_HAND_MOVEMENT_Barney

clc;
clear ;
close all;


codes_dir = 'e:\NAVEEN_Work\Cerebellum\Codes\CER_codes';
data_dir  = 'e:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS';


disp('!!!!!  ANALYSE_HAND_MOVEMENT has started running  !!!!!')
cd(data_dir)

% LOADING FILE--------------------------------------------------------
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE XML DATA FILE')
[FileName1,PathName1] = uigetfile('*.xml','File to Append');   % Open standard dialog box for retrieving files
XMLfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));

cd(PathName1);



NOME = FileName1(1:13);
TRIAL = FileName1(25);
% TRIAL = input('Enter the trial number here: ');


cd('e:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS')
mkdir(strcat(NOME));
Results_dir = strcat('e:\NAVEEN_Work\Cerebellum\Data\videos\DETAILS\',NOME);












[tracks, md] = importTrackMateTracks(XMLfile);





track_temp = cell2mat(tracks);
track_temp = (track_temp(:,1:3));

for i=2:3
  track_temp(:,i)=smooth(smooth(track_temp(:,i)));  
end


F = figure();
% plot((track_temp(:,3)),(track_temp(:,2)))
map = makeColorMap([1 0 0],[0 1 0],[0 0 1],100);
X = track_temp(:,3); Y = track_temp(:,2);
colour_line_n(X,Y,3,parula)
axis off;
colorbar
filename = strcat(NOME,'_',num2str(TRIAL));
cd(Results_dir);
print(F, '-dpdf', filename, '-r400');

clear F;





% F = figure();
% 
% X = track_temp(:,3); Y = track_temp(:,2);
% plot((track_temp(:,3)),(track_temp(:,2)),'-r','linewidth',2)
% axis off;
% filename = strcat(NOME,'_',num2str(TRIAL),'_SAMECOLOR');
% cd(Results_dir);
% print(F, '-dpdf', filename, '-r400');
% 
% clear F;








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





end