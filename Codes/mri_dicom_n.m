
%% mri_dicom_n

clc;
clear ;
close all;


codes_dir = fullfile('E:','NAVEEN_Work','Cerebellum','Codes','CER_codes_NEW','LINEAR');
% data_dir  = fullfile('E:\NAVEEN_Work\Cerebellum\Monkeys MRI\Barney\12_14_2018\data');
data_dir = fullfile('E:\NAVEEN_Work\Cerebellum\Data\MRI\Silas\9_20_17\A');
% data_dir  = fullfile('E:\NAVEEN_Work\Cerebellum\Data\MRI\Astro\11_6_19');
cd(data_dir)

%% LOAD DATA -------------------------


SIZE = 512;


cd(data_dir);
D = dir;
TEMP = extractfield(D,'name');
FILE_NAMES = (TEMP)'; FILE_NAMES = FILE_NAMES(3:end,1);
DATA = zeros(SIZE,SIZE,length(FILE_NAMES));


count = 0;
for kk=1:length(FILE_NAMES)
    str = cell2mat(FILE_NAMES(kk,:));
    if (length(str)>2)
        if (strcmp(str(1),'Z'))
            count = count+1;
            FILES(count,1:length(FILE_NAMES(kk,:) )) = FILE_NAMES(kk,:);
            try
            DATA(:,:,count) = dicomread(cell2mat(FILE_NAMES(kk,:)));
            end
        end
    end
end




  DATA_n = DATA;

% NN=50;
% DATA_n(find(DATA_n<NN))=0;


DATA_S = DATA_n;
DATA_C = permute(DATA_n,[2 3 1]);



%% VISUALIZE SINGLE SLIDE -------------------------


 i=16;
%  i=119;

F =figure();
I = (DATA_C(:,:,i));
% imagesc(imrotate(I,90),[0 25]);
imagesc(imrotate(I,90));
colormap(gray)
axis off;



filename = strcat('MRI_b_',num2str(i));
cd('E:\NAVEEN_Work\Cerebellum\Paper\1_ONE');
print(F, '-dpdf', filename, '-r1200');
clear F;



% % % 
% % % F =figure();
% % % I = (DATA_C(:,:,i));
% % % imagesc(imrotate(I,-180),[0 25]);
% % % colormap(gray)
% % % axis off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 J=122-4;

F =figure();
i=J;
subplot(2,4,1)
I = (DATA_S(:,:,i));
% imagesc(imrotate(I,90),[0 25]);
imagesc(imrotate(I,90));
colormap(gray)
axis off;
xlim([70 110])
xlim([50 90])
ylim([75 160])
ylim([60 140])
text(90,85,strcat('ML ',{' '},num2str((i-104)/2)),'fontsize',8,'color','w');

i=J+4;
subplot(2,4,2)
I = (DATA_S(:,:,i));
% imagesc(imrotate(I,90),[0 25]);
imagesc(imrotate(I,90));
colormap(gray)
axis off;
xlim([70 110])
ylim([75 160])
text(90,85,strcat('ML ',{' '},num2str((i-104)/2)),'fontsize',8,'color','w');

i=J+4;
subplot(2,4,3)
I = (DATA_S(:,:,i));
% imagesc(imrotate(I,90),[0 25]);
imagesc(imrotate(I,90));
colormap(gray)
axis off;
xlim([70 110])
ylim([75 160])
text(90,85,strcat('ML ',{' '},num2str((i-104)/2)),'fontsize',8,'color','w');

i=J+6;
subplot(2,4,4)
I = (DATA_S(:,:,i));
% imagesc(imrotate(I,90),[0 25]);
imagesc(imrotate(I,90));
colormap(gray)
axis off;
xlim([70 110])
ylim([75 160])
text(90,85,strcat('ML ',{' '},num2str((i-104)/2)),'fontsize',8,'color','w');


filename = strcat('MRI_b_FOUR2');
cd('E:\NAVEEN_Work\Cerebellum\Paper\1_ONE');
print(F, '-dpdf', filename, '-r1200');
clear F;



%%

F  = figure()
addd = 1;
START = 120+addd;
END = 124+addd;

count=0;
N = round(sqrt((END-START)));
for i=START:END
    count=count+1;
    subplot(N,N,count)
    I = (XXX(:,:,i));
    imagesc(imrotate(I,90));
    colormap(gray)
    axis off;
end



%%%%%%%%%%%%%%%%%%%%

xxx=permute(XXX,[2 3 1]);
I=xxx(:,:,125);
imagesc(I), colormap(gray)
imagesc(I'), colormap(gray)
imcontrast