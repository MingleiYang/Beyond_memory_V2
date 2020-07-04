function [Index, data_with,MeanBdisc,MeanAMdisc, P1APmean, FitP1A_All, FitP1P_All,Varlist] = ...
    plot_option_6(parameter1,parameter2,p3ser,p3_3ser,p4_3ser,p4ser,p5ser,mser,tfser,...
    Cp_eser1,Cp_eser2,Cp_eser3,Ce_pser2,Ce_pser3,P3Val,p5Val, n,...
    BIG_meanA,BIG_meanM,BIG_meanB,numstepeq,numstep1,numstep2,numstep3,numstep4,numstep5,numstep6,numstep,note)

% plot_option_6  eye disc parameter space analysis. Plots fit to eya data, bxd data and variegation.    
% Leonie Ringrose, 17.05.19


% extract data from output BIG_meanA, BIG_meanM and BIG_meanB for plots. 
%-------------------------------------------------------------------------------------
BIG_meanAM = BIG_meanA - BIG_meanM; % This converts the scale for the PRE to [-1 : +1]

eyestart = numstepeq+numstep1+numstep2+numstep3+1;

% define the last part of the time course for the first rep) 

Brep456list = BIG_meanB (eyestart:numstep, :);
AMrep456list = BIG_meanAM (eyestart:numstep, :);

Brep456sum = BIG_meanB (eyestart:numstep, :);
AMrep456sum = BIG_meanAM (eyestart:numstep, :);

x = numstep; % extract the eye disc sections and add them together

    for i = 1:n-1
   
    Brep456list = [Brep456list; (BIG_meanB ((i*x)+eyestart:(i+1)*x,:))]; % concantenates all repetitions vertically
    AMrep456list = [AMrep456list; (BIG_meanAM ((i*x)+eyestart:(i+1)*x,:))];
    
     % Add all values into the same matrix
     
    Brep456sum = Brep456sum + (BIG_meanB ((i*x)+eyestart:(i+1)*x,:));
    AMrep456sum = AMrep456sum + (BIG_meanAM ((i*x)+eyestart:(i+1)*x,:));
   
    
    end
    
    a = size (Brep456sum);
    b = a (1); % size of the time course
    c = numel(parameter1)*numel(parameter2); %81
   
    Breplist = zeros (b,n*2);
  
    for k = 1:c
        
        hic = 2*k;
        loc = hic-1;
        
        for j = 1:n-1
               
        hi = j*b;
        lo = j*b-b+1;
      
        
        Breplist (:,2*j-1:2*j) = Brep456list (lo:hi,loc:hic); %take two columns of each time block %2*k-1:2*k
        
        MeantimeB = mean (Breplist');
        STDtimeB = std (Breplist'); % calculates the standard deviation across all repetitions
        
        %numstepeye= numstep4 + numstep5 + numstep6;
        
        STDtimeBS = STDtimeB (60:180); % approximately the time points at which the fitting is performed
        Var = mean (STDtimeBS);
           
        end
        
        Varlist (k) = Var;
        
    end
    
    %take the mean
    
    MeanBdisclong = Brep456sum/n; 
    MeanAMdisclong = AMrep456sum/n;
    
    numstepeye= numstep4 + numstep5 + numstep6;
    
    % extract the last 190 rows of the data. 
   
    MeanBdisc = MeanBdisclong (numstep4-21:numstepeye,:);    %takes the last 21 windows of rep4, making a matrix with 190 rows
    MeanAMdisc = MeanAMdisclong (numstep4-21:numstepeye,:);  %takes the last 21 windows of rep4,making a matrix with 190 rows
   
  % Extract positions from mean time course (converts 190 positions to 40). 
  
    MeanBS0 = MeanBdisc (6:20,:);
    MeanAMS0 = MeanAMdisc (6:20,:);  %positions 1:15

for i = 1:9
    MeanBS1 (i,:) = MeanBdisc ((i*5+20),:);
    MeanAMS1 (i,:) = MeanAMdisc ((i*5+20),:); %positions 16:24
end

    MeanBS2 = MeanBdisc ((72:73),:);
    MeanAMS2 = MeanAMdisc ((72:73),:);  % positions 25:26

for j = 1:14                            %positions 27:40
    MeanBS3 (j,:) = MeanBdisc ((j*5+70),:);
    MeanAMS3 (j,:) = MeanAMdisc ((j*5+70),:);
end

MeanBS = [MeanBS0; MeanBS1; MeanBS2; MeanBS3];
MeanAMS = [MeanAMS0; MeanAMS1; MeanAMS2; MeanAMS3];
    

% Do the fitting to the eya PRE/TRE
%______________________________________________________________________________________________________________
% This is the raw data for eya PRE/TRE
eye_no_PRE_raw = [0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0120 0.0240  0.0480  0.1080    0.1800    0.2640    0.3960    0.4680    0.5400 0.5760  0.6000   0.5760 0.4800 0.4200  0.4680 0.4200  0.3600 0.3000  0.2640  0.2340 0.2280  0.2220 0.2160  0.2100  0.2040 0.1980  0.1920  0.1860  0.1704];
eye_with_PRE_raw = [0.01 0.01 0.01 0.01	0.01 0.01 0.01 0.01	0.01 0.01	0.01	0.01	0.02	0.04	0.09	0.15	0.28	0.41	0.54	0.67	0.78	0.8	0.78	0.5	0.4	0.39	0.3	0.2	0.16	0.14	0.12	0.105	0.09	0.075	0.06	0.045	0.03	0.02	0.015	0.01];

eye_no_PRE = (eye_no_PRE_raw)+0.1;
eye_with_PRE = (eye_with_PRE_raw)+0.1;

data_with  = eye_with_PRE'; % transpose to column

data_bxd_fit = [ 0.0978  0.1006  0.0954  0.0923  0.0907  0.0983  0.0924 0.0905  0.0852  0.0810  0.0780 0.0799   0.0870 0.0950 0.0958 0.3321 0.4351 0.5102 0.5412 0.5561 0.5760  0.5768  0.5717 0.5744 0.5579 0.5505 0.4566 0.3762 0.3801 0.3673 0.3555 0.3479 0.3616 0.3620 0.3489 0.3473 0.3627 0.3511 0.3506 0.3606]; 
% This is the best fit of bxd PRE with model 1. 

% if you only want to fit parts of the time course, then at this point,
% make a smaller version. e.g., positions 15-24 and 27-40. 

data_with_fit1 = data_with (15:24);
data_with_fit2 = data_with (29:40);
data_with_fit = [data_with_fit1;data_with_fit2];

MeanBS_fit1 = MeanBS (15:24,:);
MeanBS_fit2 = MeanBS (29:40,:);
MeanBS_fit = [MeanBS_fit1;MeanBS_fit2];% (40,81*2)

   
    z = numel(parameter1)*numel(parameter2)*2;
    
    for i = 1:z
        controlMatrix (:,i) = data_with_fit; % makes a matrix will all columns identical (40,81)                                                                                                                                                             
        controlMatrix_bxd (:,i) = data_bxd_fit';
    
    end
    
    Bdisc_diff_all = MeanBS_fit - controlMatrix; % subtracts the data values from control matrix from each time coruse
    Bdisc_diff_all2 = (Bdisc_diff_all).^2;      % square of each position.    
    MeanBdisc_diff2 = mean (Bdisc_diff_all2);   % mean of each column.
  
    bxd_diff_all = MeanBS - controlMatrix_bxd;
    bxd_diff_all2 = (bxd_diff_all).^2;
    Mean_bxd_diff2 = mean (bxd_diff_all2);
    
    sum = numel (MeanBdisc_diff2)/2;
    
    P1A_All = MeanBdisc_diff2 (1); %extracts the 1st position for P1 input for the "anterior" 
    P1P_All = MeanBdisc_diff2 (2); %extracts the 2nd position for P1 input for the "posterior" (in thsi case p1 is the same in anterior and posterior, they are treated as replicates.
    
    bxdA_All = Mean_bxd_diff2 (1); %extracts the 1st position for P1 input for the "anterior" 
    bxdP_All = Mean_bxd_diff2 (2); %extracts the 2nd position for P1 input for the "posterior" (in thsi case p1 is the same in anterior and posterior, they are treated as replicates.
    
    
    for i = 1:sum-1
        
        P1A_All = [P1A_All,MeanBdisc_diff2(2*i+1)]; % every odd numbered position
        P1P_All = [P1P_All,MeanBdisc_diff2(2*i+2)]; % every even numbered position
        
        bxdA_All = [bxdA_All,Mean_bxd_diff2(2*i+1)]; % every odd numbered position
        bxdP_All = [bxdP_All,Mean_bxd_diff2(2*i+2)]; % every even numbered position
    end 
    
    P1AP_All = [P1A_All;P1P_All]; % Puts them in a matrix of 2 rows
    P1APmean = mean (P1AP_All); % takes the mean of both rows. This is the vector of fit values for each parameter combination (1,81)
    [Min,Index] = min (P1APmean); %extracts the value and position of the best fit. Then this is used for the plot. 
    % can also extract the parameter values for that: p3ser (min), p5ser
    % (min).
    
    bestp3 = P3Val (Index);
    bestp5 = p5Val (Index);
    bestp3p5 = [bestp3;bestp5];
    
    bxdAP_All = [bxdA_All;bxdP_All]; % Puts them in a matrix of 2 rows
    bxdAPmean = mean (bxdAP_All); % takes the mean of both rows. This is the vector of fit values for each parameter combination (1,81)
    [Minb,Indexb] = min (bxdAPmean); %extracts the value and position of the best fit. Then this is used for the plot. 
    % can also extract the parameter values for that: p3ser (min), p5ser
    % (min).
    
    bestp3bxd = P3Val (Indexb);
    bestp5bxd = p5Val (Indexb);
    bestp3p5bxd = [bestp3bxd;bestp5bxd];
    
   % Extract the data for the plot axes titles from the input vectors. 
%-------------------------------------------------------------------------------------------------
% inputs: parameter1, parameter2, n, V, C, BIG_meanB, BIG_meanA

if parameter1 == p3ser
    parameter1_label = 'Feedback (p3,p4)';
    X2 = 'Feedback (p3,p4)';
    name1 = '(p3,p4)';
    
elseif parameter1 == p4ser
    parameter1_label = 'parameter 1: p 4 (bias to A)';
    X2 = 'parameter 1: p4';
    name1 = 'p4';
    
elseif parameter1 == p3_3ser
    parameter1_label = 'parameter 1: p3 in PcG mutant';
    X2 = 'parameter 1: p3_3';
    name1 = 'p3_3';
    
elseif parameter1 == p4_3ser
    parameter1_label = 'parameter 1: p4 in TrxG mutant';
    X2 = 'parameter 1: p4_3';
    name1 = 'p4_3';
    
elseif parameter1 == p5ser
    parameter1_label = 'parameter 1: p 5 (noise)';
    X2 = 'parameter 1: p5';
    name1 = 'p5';

elseif parameter1 == mser
    parameter1_label = 'parameter 1: m (PRE/TRE size)';
    X2 = 'parameter 1: m';
    name1 = 'm';
    
elseif parameter1 == tfser
    parameter1_label = 'parameter 1: tf (promoter binding sites)';
    X2 = 'parameter 1: tf; ';
    name1 = 'tf';
    
elseif parameter1 == Cp_eser1
    parameter1_label = 'parameter 1: Cp_e1 (PRE/TRE to promoter)';
    X2 = 'parameter 1: Cp_e1; ';
    name1 = 'Cp_e1';
    
elseif parameter1 == Cp_eser2
    parameter1_label = 'parameter 1: Cp_e2 (PRE/TRE to promoter)';
    X2 = 'parameter 1: Cp_e2; ';
    name1 = 'Cp_e2';
    
elseif parameter1 == Cp_eser3
    parameter1_label = 'parameter 1: Cp_e3 (PRE/TRE to promoter)';
    X2 = 'parameter 1: Cp_e3; ';
    name1 = 'Cp_e3';
    
elseif parameter1 == Ce_pser2
    parameter1_label = 'parameter 1: Ce_p2 (promoter to PRE/TRE)';
    X2 = '(parameter 1: Ce_p2; ';
    name1 = 'Ce_p2';

elseif parameter1 == Ce_pser3
    parameter1_label = 'parameter 1: Ce_p3 (promoter to PRE/TRE)';
    X2 = '(parameter 1: Ce_p3; ';
    name1 = 'Ce_p3';
    
else
    error 'parameter 1 definition'
end


if parameter2 == p3ser
    parameter2_label = 'parameter 2: p 3 (bias to M)';
    X1 = 'parameter 2: p3';
    name2 = 'p3';
    
elseif parameter2 == p4ser
    parameter2_label = 'parameter 2: p 4 (bias to A)';
    X1 = 'parameter 2: p4';
    name2 = 'p4';
    
elseif parameter2 == p3_3ser
    parameter2_label = 'parameter 2: p 3 in PcG mutant';
    X1 = 'parameter 2: p3_3';
    name2 = 'p3_3';
    
elseif parameter2 == p4_3ser
    parameter2_label = 'parameter 2: p 4 in TrxG mutant';
    X1 = 'parameter 2: p4_3';
    name2 = 'p4_3';
    
elseif parameter2 == p5ser
    parameter2_label = 'Feedback - independent transitions (p5)';
    X1 = 'parameter 2: p5';
    name2 = 'p5';

elseif parameter2 == mser
    parameter2_label = 'parameter 2: m (PRE/TRE size)';
    X1 = 'parameter 2: m';
    name2 = 'm';

elseif parameter2 == tfser
    parameter2_label = 'parameter 2: tf (promoter binding sites)';
    X1 = 'parameter 2: tf; ';
    name2 = 'tf';
    
elseif parameter2 == Cp_eser1
    parameter2_label = 'parameter 2: Cp_e1 (PRE/TRE to promoter)';
    X1 = 'parameter 2: Cp_e1; ';
    name2 = 'Cp_e1';
    
elseif parameter2 == Cp_eser2
    parameter2_label = 'parameter 2: Cp_e2 (PRE/TRE to promoter)';
    X1 = 'parameter 2: Cp_e2; ';
    name2 = 'Cp_e2';
    
elseif parameter2 == Cp_eser3
    parameter2_label = 'parameter 2: Cp_e3 (PRE/TRE to promoter)';
    X1 = 'parameter 2: Cp_e3; ';
    name2 = 'Cp_e3';
    
elseif parameter2 == Ce_pser2
    parameter2_label = 'parameter 2: Ce_p2 (promoter to PRE/TRE)';
    X1 = '(parameter 2: Ce_p2; ';
    name2 = 'Ce_p2';

elseif parameter2 == Ce_pser3
    parameter2_label = 'parameter 2: Ce_p3 (promoter to PRE/TRE)';
    X1 = '(parameter 2: Ce_p3; ';
    name2 = 'Ce_p3';
    
else
    error 'parameter 2 definition'
end

% generate the lists of parameter values for the plot axes
%-----------------------------------------------------------------------------------------

num1 = numel (parameter1); %size of the block of parameter1 values. (will be rows)
num2 = numel (parameter2); %size of the block of parameter2 values. (will be columns)

parameter1_list = name1;
parameter2_list = name2;


for Str = 1:num1
 
    par1 = mat2str(parameter1(Str));
    parameter1_list = char(parameter1_list,par1);
    
end


for Str = 1:num2 
 
    par2 = mat2str(parameter2(Str));
    parameter2_list = char(parameter2_list,par2);
    
end

X3= 'parameter 1: ';
X4= ', ';

for par1 = 1:num1+1
    
    X3 = [X3,X4,(parameter1_list (par1,:))];
    
end


Xlabel2 = {X1;X3};

% use note and Xlabel2 in plots!    
        
FitP1A_All = reshape (P1A_All, num1,num2);
FitP1P_All = reshape (P1P_All, num1,num2);
FitP1APmean = reshape (P1APmean, num1,num2);

Var = reshape (Varlist,num1,num2);

Fit_bxd = reshape (bxdAPmean, num1,num2);

lastcol = zeros (num1,1); 
lastrow = zeros (1,num2+1);    


FitP1A_All0 =  [FitP1A_All,lastcol];
FitP1A_All00 =  [FitP1A_All0;lastrow];

  
FitP1P_All0 =  [FitP1P_All,lastcol];
FitP1P_All00 =  [FitP1P_All0;lastrow];

lastcolAP = FitP1APmean(:,num1);
FitP1AP0 =  [FitP1APmean,lastcolAP];
lastrowAP = FitP1AP0(num2,:);
FitP1AP00 =  [FitP1AP0;lastrowAP];

lastcolVar = Var(:,num1);
Var0 =  [Var,lastcolVar];
lastrowVar = Var0(num2,:);
Var00 =  [Var0;lastrowVar];

lastcolbxd = Fit_bxd(:,num1);
bxd0 =  [Fit_bxd,lastcolbxd];
lastrowbxd = bxd0(num2,:);
bxd00 =  [bxd0;lastrowbxd];

figure 
colormap (flipud (parula));

 title1 = {'Plot option 6, Eye disc parameter space';note;'';'fit to eya data'};
 
 subplot1 = subplot(3,3,1);
    pcolor (FitP1AP00); caxis ([0.005,0.011]);shading interp;
    title (title1);
    set(gca,'XTick',1:1:num2+1)
    set(gca,'XTickLabel',(parameter2_list));
    set(gca,'YTick',1:1:num1+1)
    set(gca,'YTickLabel',(parameter1_list));
    xlabel(parameter2_label);
    ylabel(parameter1_label);
    
 subplot2 = subplot(3,3,4);
    plot (data_with);
    ylim ([0 1]);    
    hold all; plot (MeanBS(:,Index*2));
    xlabel('position in disc (A-P)'); 
    ylabel('Mean promoter (N_B/s)');
    title ('best fit eya');
    legend('data ','model');

 subplot3 = subplot(3,3,7);
    bar (bestp3p5); 
    ylim ([0 0.8]);
    set(gca,'XTickLabel',{'best (p3,p4)','best p5'});
    title ('values of (p3,p4) and p5 for best eya fit');
 
 subplot4 = subplot(3,3,2);
    pcolor (bxd00); caxis ([0.004,0.01]);shading interp;
    title ('fit to bxd data');
    set(gca,'XTick',1:1:num2+1)
    set(gca,'XTickLabel',(parameter2_list));
    set(gca,'YTick',1:1:num1+1)
    set(gca,'YTickLabel',(parameter1_list));
    xlabel(parameter2_label);
    ylabel(parameter1_label);
 
 subplot5 = subplot(3,3,5);
    plot (data_bxd_fit);
    ylim ([0 1]);    
    hold all; plot (MeanBS(:,Indexb*2));
    xlabel('position in disc (A-P)'); 
    ylabel('Mean promoter (N_B/s)');
    title ('best fit bxd');
    legend('data','model');

 subplot6 = subplot(3,3,8);
    bar (bestp3p5bxd); %bestp3p5% for all fits, plot P1APmean   ylim ([0 0.05]);
    ylim ([0 0.8]);
    set(gca,'XTickLabel',{'best (p3,p4)','best p5'});
    title ('values of (p3,p4) and p5 for best bxd fit');
 
 subplot7 = subplot(3,3,3);
    pcolor (Var00); caxis ([0.2,0.35]);shading interp;
    title ('variegation test');
    set(gca,'XTick',1:1:num2+1)
    set(gca,'XTickLabel',(parameter2_list));
    set(gca,'YTick',1:1:num1+1)
    set(gca,'YTickLabel',(parameter1_list));
    xlabel(parameter2_label);
    ylabel(parameter1_label); 
 
bestp3 
bestp5

bestp3bxd
bestp5bxd
end

