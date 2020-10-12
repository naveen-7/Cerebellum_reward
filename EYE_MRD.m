
%% ANALYSE_EYE with MRD
%% Written by naveen at JLG on 3/14/19


clc;
clear all;
close all;


% Setup directories--------------------------------------------------------

codes_dir = fullfile('e:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
data_dir  = fullfile('e:','NAVEEN_Work','Cerebellum','Data','MERGED_CELLS');

cd(data_dir)

% LOADING FILE------------------------------------------------------- %% 
disp('*************************************************')
disp('*************************************************')
disp('LOAD THE CELL DATA FILE')
[FileName1,PathName1] = uigetfile('*.mat','File to Append');   % Open standard dialog box for retrieving files
DATAfile = cat(2,PathName1,FileName1);
disp(strcat('!!!!!','File you entered is :',FileName1,' !!!!!'));
load(DATAfile);
cd(PathName1);


%% CASE

F = figure();

TEMP_infos = Infos; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Infos = Infos; CASE = '1-N'; COL = 3; % 1-N

% Infos(2:end,10) = Infos(1:end-1,10); CASE = '1+N'; COL = 4; % 1+N
% Infos(1:end-1,10) = Infos(2:end,10); CASE = '2-N'; COL = 2; % 2-N

  Infos(3:end,10) = Infos(1:end-2,10); CASE = '2+N'; COL = 5; % 2+N
% Infos(1:end-2,10) = Infos(3:end,10); CASE = '3-N'; COL = 1; % 3-N

NN = 70;
INFOS_CHANGE = Infos(CHANGE:CHANGE+NN,:);
EYEX_CHANGE = EYE.EYE_Xn(CHANGE:CHANGE+NN,:);
EYEY_CHANGE = EYE.EYE_Yn(CHANGE:CHANGE+NN,:);

IND = find(~cellfun(@isempty,EYEX_CHANGE));
INFOS_CHANGE = INFOS_CHANGE(IND,:);
EYEX_CHANGE  = EYEX_CHANGE(IND,:);
EYEY_CHANGE = EYEY_CHANGE(IND,:);




% total signal = @M(N) (200:700) + @S(N+1) (-400:200) + @M(N+1) (-300: 200)

%%%%%% part1: @M(N) (200:700)
clear PART1X PART1Y PART1 TIME
START = 200; END = 700; TIME = START:END;

PART1X = NaN(size(INFOS_CHANGE,1)-1,length(TIME));
PART1Y = NaN(size(INFOS_CHANGE,1)-1,length(TIME));

for i=1:size(INFOS_CHANGE,1)-1
    PART1X(i,:) = EYEX_CHANGE{i,1}(INFOS_CHANGE(i,11)+START:INFOS_CHANGE(i,11)+END);
    PART1Y(i,:) = EYEY_CHANGE{i,1}(INFOS_CHANGE(i,11)+START:INFOS_CHANGE(i,11)+END);
end
PART1 = [PART1X PART1Y];



%%%%%% part2: @S(N+1) (-300:200)
clear PART2X PART2Y PART2 TIME
START = -300; END = 200; TIME = START:END;

PART2X = NaN(size(INFOS_CHANGE,1)-1,length(TIME));
PART2Y = NaN(size(INFOS_CHANGE,1)-1,length(TIME));

for i=1:size(INFOS_CHANGE,1)-1
    PART2X(i,:) = EYEX_CHANGE{i+1,1}(INFOS_CHANGE(i+1,4)+START:INFOS_CHANGE(i+1,4)+END);
    PART2Y(i,:) = EYEY_CHANGE{i+1,1}(INFOS_CHANGE(i+1,4)+START:INFOS_CHANGE(i+1,4)+END);
end
PART2 = [PART2X PART2Y];



%%%%%% part3: @M(N+1) (-300: 200)
clear PART3X PART3Y PART3 TIME
START = -300; END = 200; TIME = START:END;

PART3X = NaN(size(INFOS_CHANGE,1)-1,length(TIME));
PART3Y = NaN(size(INFOS_CHANGE,1)-1,length(TIME));

for i=1:size(INFOS_CHANGE,1)-1
    if INFOS_CHANGE(i+1,11)>300
        PART3X(i,:) = EYEX_CHANGE{i+1,1}(INFOS_CHANGE(i+1,11)+START:INFOS_CHANGE(i+1,11)+END);
        PART3Y(i,:) = EYEY_CHANGE{i+1,1}(INFOS_CHANGE(i+1,11)+START:INFOS_CHANGE(i+1,11)+END);
    end
end
PART3 = [PART3X PART3Y];


ALLPARTS = [PART1 PART2 PART3];
corrind = find(INFOS_CHANGE(1:size(ALLPARTS,1)-1,10)==1);
wrngind = find(INFOS_CHANGE(1:size(ALLPARTS,1)-1,10)==2);

ALLPARTS_C = ALLPARTS(corrind,:);
ALLPARTS_W = ALLPARTS(wrngind,:);


F = figure();
subplot(2,2,1)
[p,~,MRD_a, MRD_w] = MRD_n(ALLPARTS_C,ALLPARTS_W,250,0.8,1);



subplot(2,2,2)
hold on;
errorline_n(1:length(ALLPARTS_C),nanmean(ALLPARTS_C),nanstd(ALLPARTS_C)/sqrt(size(ALLPARTS_C,1)),'1',[0 0 1]);
errorline_n(1:length(ALLPARTS_W),nanmean(ALLPARTS_W),nanstd(ALLPARTS_W)/sqrt(size(ALLPARTS_W,1)),'1',[1 0 0]);

for i=1:size(ALLPARTS_C,2)
PP(i) = ttest_NN(ALLPARTS_C(:,i),ALLPARTS_W(:,i));
end

pind = find(PP<0.05);
plot(pind,20,'.k')


cd(Results_dir)
filename = 'EYE_MRD';
print(F, '-dpdf', filename, '-r400')

