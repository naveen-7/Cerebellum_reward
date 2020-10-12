% function to plot a scatter plot and a overlayed histogram showing x and y
% distributions

% created by naveen at cumc on 10/23/17

function scatter_hist_n(X,Y,MS,Colour,TURN)


% X*             : N x 1   : array of x values
% Y*             : N x 1   : array of x values
% MS             : 1 x 1   : start time of the PSTH
% Colour         : 1 x 3   : RGB value of raster   [Default: black]
% TURN           : 1 x 1   : 1 to turn the y plot, 0 not to turn [Default: 1]


if nargin<2
    error('Incomplete input to the function Raster_n');
elseif nargin==3
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = MS;
    Colour      = [0 0 0];
    TURN=1;
elseif nargin==4
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = MS;
    varargin{4} = Colour;
    TURN=1;
elseif nargin==5
    varargin{1} = X;
    varargin{2} = Y;
    varargin{3} = MS;
    varargin{4} = Colour;
    varargin{5} = TURN;
else
    error('Too many inputs to the function scatter_hist_n');
end


BIN = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S1 = subplot(10,10,[12 13 14 15 22 23 24 25 32 33 34 35 42 43 44 45]);
S1 = subplot(10,10,[12 13 14 15 22 23 24 25 32 33 34 35 42 43 44 45],'box','off','tickdir','out','Linewidth',0.3,'FontSize',7);

S1Pos = get(S1,'position');
hold on;
plot(X, Y,'o','markersize',MS,'color',Colour)
plot(nanmean(X), nanmean(Y),'O','markersize',7,'color','k','markerfacecolor',Colour)

MIN = round(nanmin([X;Y]));
MAX = round(nanmax([X;Y]));
xlim([MIN MAX]); ylim([MIN MAX])
plot([MIN MAX],[MIN MAX],'--','color',[0.7 0.7 0.7])
set(gca,'fontsize',8,'linewidth',0.3)





S2 = subplot(10,10,[2 3 4 5],'box','off','tickdir','out','Linewidth',0.3,'FontSize',7);

S2Pos = get(S2,'position');
hold on;
h = histogram(X,BIN);
h.FaceColor = Colour;
h.EdgeColor = 'w';
box off;
xlim([MIN MAX]);
set(gca,'fontsize',8,'linewidth',0.3)
set(gca,'xtick',[])
set(gca,'xticklabel',[])


S3 = subplot(10,10,[11 21 31 41],'box','off','tickdir','out','Linewidth',0.3,'FontSize',7);

S3Pos = get(S3,'position');
hold on;
h = histogram(Y,BIN);
h.FaceColor = Colour;
h.EdgeColor = 'w';
box off;
xlim([MIN MAX]);
set(gca,'fontsize',8,'linewidth',0.3)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
if TURN==1
    camroll(90)
end



set(S1,'position',S1Pos);
set(S2,'position',S2Pos);
set(S3,'position',S3Pos);


end