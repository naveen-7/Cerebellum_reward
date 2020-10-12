%% Function to correlate trails with themselves, for state change analysis
% Created by Naveen on 3/29/17 at cumc


function BIG_CORR_n(Matrix,LC)

% Matrix is the matrix of PSTH values in that epoch
% LC is a 2 x n matrix whose first row is x values and second row is %corr

if nargin<1
    error('Incomplete input to the function BIG_CORR_n');
elseif nargin==1
    varargin{1} = Matrix;
    isLC        = 0;
elseif nargin==2
    varargin{1} = Matrix;
    varargin{2} = LC;
    isLC        = 1;
else
    error('Too many inputs to the function BIG_CORR_n');
end





%% 

Coeff = NaN(size(Matrix,1),size(Matrix,1));
COEFF = NaN(size(Matrix,1),size(Matrix,1));
PVal = NaN(size(Matrix,1),size(Matrix,1));
Pval = NaN(size(Matrix,1),size(Matrix,1));
PVAL = NaN(size(Matrix,1),size(Matrix,1));
for i = 1:size(Matrix,1)
    TEMP1 = Matrix(i,:);
    for j = 1:size(Matrix,1)
        TEMP2 = Matrix(j,:);
         [Coeff(i,j),PVal(i,j)] = corr(TEMP1',TEMP2');
       if PVal(i,j)<=0.05
           COEFF(i,j)=Coeff(i,j);
       end
         Pval(i,j) = ranksum(TEMP1',TEMP2');
        if Pval(i,j)<0.05
            PVAL(i,j)=0;
        end
        if Pval(i,j)>=0.05
            PVAL(i,j)=1;
        end
        
    end
end


PVAL(find(Pval<0.05))=1;
PVAL(find(Pval>=0.05))=0;

% Coeff(find(isnan(Coeff)))=1;




PLOTTHIS = COEFF;
figure
imagesc(PLOTTHIS); 
colormap(jet)
% colormap(gray)
colorbar
% caxis([0 1]);






end