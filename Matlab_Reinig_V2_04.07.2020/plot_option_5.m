function [Timeno, Timewith, TimewithS,Discno, Discwith, AMDiscno, AMDiscwith, STDno, STDwith] = ...
    plot_option_5(BIG_meanA,BIG_meanM,BIG_meanB,n,...
    numstepeq,numstep1,numstep2,numstep3,numstep4,numstep5,numstep6,numstep, note, PRE_ID)

% plot_option_5  eye disc. Extracts the data for third instar from whole time course. Compiles eye disc images and plots average profiles. 
% Leonie Ringrose, 17.05.19


% extract data from output BIG_meanA, BIG_meanM and BIG_meanB for plots. 
%------------------------------------------------------------------------------------------------------------------

BIG_meanAM = BIG_meanA - BIG_meanM; % This converts the scale for the PRE to [-1 : +1]

eyestart = numstepeq+numstep1+numstep2+numstep3+1;

% define the last part of the time course for the first rep) 

Brep456 = BIG_meanB (eyestart:numstep, :);
AMrep456 = BIG_meanAM (eyestart:numstep, :);


x = numstep; % extract the eye disc sections and join them vertically


    for i = 1:n-1
   
    Brep456 = [Brep456; (BIG_meanB ((i*x)+eyestart:(i+1)*x,:))];
    AMrep456 = [AMrep456; (BIG_meanAM ((i*x)+eyestart:(i+1)*x,:))];
    
    end
    
    % Each matrix has 4 columns, [left, right], without, with coupling. n repetitions are concatenated vertically
    % Reshape each column into a matrix of n columns, with one row for each window. (numstepeye)
    
    numstepeye= numstep4 + numstep5 + numstep6;
    
    Brep456noL = reshape (Brep456 (:,1), numstepeye,n); 
    Brep456noR = reshape (Brep456 (:,2), numstepeye,n); 
    
    Brep456withL = reshape (Brep456 (:,3), numstepeye,n); 
    Brep456withR = reshape (Brep456 (:,4), numstepeye,n); 
   
    % join left and right. 
    
    
    DiscnoALL = [Brep456noL,Brep456noR];
    DiscwithALL = [Brep456withL,Brep456withR];
    
    Discno = DiscnoALL (numstep4-21:numstepeye,:);      %takes the last 21 windows of rep4
    Discwith = DiscwithALL (numstep4-21:numstepeye,:);  %takes the last 21 windows of rep4
  
    % calculate means for time course fig. Currently calculates for whole
    % of rep4. 
    
    Timeno = mean (Discno');
    STDno = std (Discno');
    
        
    Timewith = mean (Discwith');
    STDwith = std (Discwith');
    
    % repeat for the PRE/TRE
    
    AMrep456noL = reshape (AMrep456 (:,1), numstepeye,n); 
    AMrep456noR = reshape (AMrep456 (:,2), numstepeye,n); 
    
    AMrep456withL = reshape (AMrep456 (:,3), numstepeye,n); 
    AMrep456withR = reshape (AMrep456 (:,4), numstepeye,n); 
    
    % join left and right. 
    
    AMDiscnoALL = [AMrep456noL,AMrep456noR];
    AMDiscwithALL = [AMrep456withL,AMrep456withR];
  
    AMDiscno = AMDiscnoALL (numstep4-21:numstepeye,:); %takes the last 21 windows of rep4
    AMDiscwith = AMDiscwithALL (numstep4-21:numstepeye,:); %takes the last 21 windows of rep4 
  
    % calculate means for time course fig. 
    
    AMTimeno = mean (AMDiscno');
    AMSTDno = std (AMDiscno');
    
        
    AMTimewith = mean (AMDiscwith');
    AMSTDwith = std (AMDiscwith');
    
    %Calculate mean promoter state for the end of the time course (last 20
    %windows, = 200 mins). This is a proxy for adult eye colour
    %---------------------------------------------------------------------------------------------------------------------
    
    Dn = size (Discno);
    Dnlast = Dn (1);
    Dw = size (Discwith);
    Dwlast = Dw (1);
    
    A1 = Discno (Dnlast-20:Dnlast,:);
    A2 = Discwith (Dwlast-20:Dwlast,:);
    
    A1m = mean (A1);
    A2m = mean (A2);
    
    meanA1m = mean (A1m);
    meanA2m = mean (A2m);
    
    
    stdA1m = std (A1m);
    stdA2m = std (A2m);
    
    Enhancermeans = [meanA1m,meanA2m];
    Enhancerstds = [stdA1m,stdA2m];
   
    PRETimewith = (AMTimewith/2) + 0.5;
    PRETimeno = (AMTimeno/2) + 0.5;
    

% calculate the fitted time course for the promoter and the PRE/TRE (this fits the raw data to the experimental data by assuming different cell sizes in zones 2 and 3)
TimenoS0 = Timeno (6:20);
TimewithS0 = Timewith (6:20);  %positions 1:15

PRETimenoS0 = PRETimeno (6:20);
PRETimewithS0 = PRETimewith (6:20);  %positions 1:15

for i = 1:9
    TimenoS1 (i) = Timeno (i*5+20);
    TimewithS1 (i) = Timewith (i*5+20); %positions 16:24
    
    PRETimenoS1 (i) = PRETimeno (i*5+20);
    PRETimewithS1 (i) = PRETimewith (i*5+20); %positions 16:24
end

TimenoS2 = Timeno (72:73);
TimewithS2 = Timewith (72:73);  % positions 25:26

PRETimenoS2 = PRETimeno (72:73);
PRETimewithS2 = PRETimewith (72:73);  % positions 25:26

for j = 1:14                            %positions 27:40
    TimenoS3 (j) = Timeno (j*5+70);
    TimewithS3 (j) = Timewith (j*5+70);
    
    PRETimenoS3 (j) = PRETimeno (j*5+70);
    PRETimewithS3 (j) = PRETimewith (j*5+70);
end

TimenoS = [TimenoS0, TimenoS1, TimenoS2, TimenoS3];
TimewithS = [TimewithS0, TimewithS1, TimewithS2, TimewithS3];

PRETimenoS = [PRETimenoS0, PRETimenoS1, PRETimenoS2, PRETimenoS3];
PRETimewithS = [PRETimewithS0, PRETimewithS1, PRETimewithS2, PRETimewithS3];

% Extract disc image (the first 140 rows of the Disc matrices).

DiscnoS = Discno (1:140,:);
DiscwithS = Discwith (1:140,:);
AMDiscnoS = AMDiscno (1:140,:);
AMDiscwithS = AMDiscwith (1:140,:);
    
% Assemble the disc from individual "cells" selected from different
% simulations at random for each position along the time course.

a = size (DiscnoS); % rows = time points, columns = repetitions (n                                                                                                                                                                                                                                                                                             *2).
b = a (1); 
c = a (2);
cell = 10; % this is the number of windows in the time course to join consecutely for each "cell"
z= int8(b/cell); %


testnoS = zeros (z,c);
testwithS = zeros (z,c);

AMtestnoS = zeros (z,c);
AMtestwithS = zeros (z,c);



for i = 1:z;  
    for j = 1:c; 
        
        r = randi (c);
        lo = cell*i-(cell-1);
        hi = cell*i; 
        
        testnoS (lo:hi,j)= DiscnoS (lo:hi,r); % takes 3 windows of the time course from a randomly chosen repetition.
        testwithS (lo:hi,j)= DiscwithS (lo:hi,r);
        
        AMtestnoS (lo:hi,j)= AMDiscnoS (lo:hi,r); 
        AMtestwithS (lo:hi,j)= AMDiscwithS (lo:hi,r);
    end
end




  

figure
colormap (bone) % This plots the data assembled from selected AP positions (to fit the data) and randomly selected individual simulations (to simluate individual cells) 

if PRE_ID  == ('e')
    title1 = {'Plot option 5, eye disc, eya PRE/TRE';note;'';'Promoter'};
elseif PRE_ID  == ('b')
    title1 = {'Plot option 5, eye disc, bxd PRE/TRE';note;'';'Promoter'};
else error 'please enter PRE_ID as either e or b'
end


subplot1 = subplot(3,2,1);
    plot (TimenoS (11:40));
    ylim ([0 1]);
    xlim ([1 30]);
    hold on; 
    plot (TimewithS (11:40));
    title (title1);
    set(gca,'XTickLabel',{'A',' ','MF',' ',' ','P'});
    ylabel('Mean promoter (N_B/s)');  
    legend('without coupling','with  coupling');

subplot2 = subplot (3,2,2);
    plot (PRETimenoS(11:40));   
    ylim ([0 1]);
    xlim ([1 30]);
    hold on; 
    plot (PRETimewithS(11:40)); 
    title ('PRE/TRE');
    set(gca,'XTickLabel',{'A',' ','MF',' ',' ','P'});
    ylabel('Mean PRE/TRE (N_A/S)');
    legend('without coupling','with  coupling');
  
subplot3 = subplot (3,2,3);
  pcolor (testnoS(12:120,1:20)'); caxis ([0,1]);shading flat;
  set(gca,'XTickLabel',{'A',' ','MF',' ','P'});
  ylabel('cell row');
  title ('Eye disc promoter without coupling');
 
subplot4 = subplot (3,2,4);
  pcolor (AMtestnoS(12:120,1:20)'); caxis ([-1,1]);shading flat; 
  set(gca,'XTickLabel',{'A',' ','MF',' ','P'});
  ylabel('cell row');
  title ('Eye disc PRE/TRE without coupling');

subplot5 = subplot (3,2,5);
  pcolor (testwithS(12:120,1:20)'); caxis ([0,1]);shading flat;
  set(gca,'XTickLabel',{'A',' ','MF',' ','P'});
  ylabel('cell row');
  title ('Eye disc promoter with coupling');
  
 subplot6 = subplot (3,2,6);
  pcolor (AMtestwithS(12:120,1:20)'); caxis ([-1,1]);shading flat; 
  set(gca,'XTickLabel',{'A',' ','MF',' ','P'});
  ylabel('cell row');
  title ('Eye disc PRE/TRE with coupling');
  
end

