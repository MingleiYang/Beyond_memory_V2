function [One_diffBm2,meanLRBm,LRmeanBm0,singlevalBm] = ...
    plot_option_2(parameter1,parameter2,p3ser,p3_3ser,p4_3ser,p4ser,p5ser,mser,tfser,...
    Cp_eser1,Cp_eser2,Cp_eser3,Ce_pser2,Ce_pser3,n,V,C,BIG_meanB,BIG_meanA,numstep,numstep3,note)

% plot_option_2  parameter space for memory of silencing in embryo, used in Figure S1 and S2 of Reinig et al.
% Leonie Ringrose, 17.05.19

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

% extract data from output BIG_meanA and BIG_meanB for plots. 
%------------------------------------------------------------------------------------------------------------------
    
for cell = 1:n
    
    embryoB1(cell,:) = BIG_meanB (cell*numstep-numstep3,:);%extract the last window of the initiation phase
    embryoB2(cell,:) = BIG_meanB (cell*numstep,:);%extract the last window of the maintenance phase

    embryoA1(cell,:) = BIG_meanA (cell*numstep-numstep3,:);%extract the last window of the initation phase
    embryoA2(cell,:) = BIG_meanA (cell*numstep,:);%extract the last window of the mantenance phase

end

controlBe = embryoB1(:,1:2); %extracts the first two columns of the promoter without coupling for initiation phase, these are the control values
controlB1 = controlBe; 


for x = 1:(V*C/2)-1 
    controlB1= [controlB1,controlBe]; 
end  

% control B1 is the promoter in early embryo without coupling
% calculate differences between experiment and control within each
% simulation.

diffBe = embryoB1 - controlB1;   %calculates the difference between promoter at end of initiation phase and control.
diffBm = embryoB2 - controlB1;   %calculates the difference between promoter at end of maintenance phase and control.  
diffAe = embryoA1 - controlB1;   %calculates the difference between PRE/TRE at end of initiation phase and control.
diffAm = embryoA2 - controlB1;   %calculates the difference between PRE/TRE at end of maintenance phase and control.  


diffBe2 = diffBe.^2;
diffBm2 = diffBm.^2;
diffAe2 = diffAe.^2;
diffAm2 = diffAm.^2;

meandiffBe = mean(diffBe2);
meandiffBm = mean(diffBm2);
meandiffAe = mean(diffAe2); 
meandiffAm = mean(diffAm2);

One_diffBm2 = 1 - diffBm2; % values in diffBm2 are between 0 and 1.
meanOne_diffBm2 = mean (One_diffBm2); 
stdOne_diffBm2 = std (One_diffBm2);


% calculate the mean and std of the promoter values in the last window of
% the maintenance phase, for Figure 3.
meanembryoB2 = mean (embryoB2);
stdembryoB2 = std (embryoB2);


for i = 1: V*C/2
    
    % extract the left and rigt values of the fit.
   
    meanLBe (i) = meandiffBe (2*i-1); % extract the left values for the promoter
    meanRBe (i) = meandiffBe (2*i); % extract the right values for the promoter 
    meanLRBe (i) = (meandiffBe (2*i-1) + meandiffBe (2*i) )/2; % calculate the mean of each pair of left and right columns (promoter)
    
    meanLBm (i) = meandiffBm (2*i-1); % extract the left values for the promoter
    meanRBm (i) = meandiffBm (2*i); % extract the right values for the promoter 
    meanLRBm (i) = (meandiffBm (2*i-1) + meandiffBm (2*i) )/2; % calculate the mean of each pair of left and right columns (promoter)
    
    meanLAe (i) = meandiffAe (2*i-1); % extract the left values for the PRE/TRE
    meanRAe (i) = meandiffAe (2*i); % extract the right values for the PRE/TRE 
    meanLRAe (i) = (meandiffAe (2*i-1) + meandiffAe (2*i) )/2; % calculate the mean of each pair of left and right columns (PRE/TRE)
    
    meanLAm (i) = meandiffAm (2*i-1); % extract the left values for the PRE/TRE
    meanRAm (i) = meandiffAm (2*i); % extract the right values for the PRE/TRE 
    meanLRAm (i) = (meandiffAm (2*i-1) + meandiffAm (2*i) )/2; % calculate the mean of each pair of left and right columns (PRE/TRE)
    
    mean_OneLRBm (i) = (meanOne_diffBm2 (2*i-1) + meanOne_diffBm2 (2*i) )/2; % calculate the mean of each pair of left and right columns (promoter)
    std_OneLRBm (i) = (stdOne_diffBm2 (2*i-1) + stdOne_diffBm2 (2*i) )/2; % calculate the mean of each pair of left and right columns (PRE/TRE)

    singlevalBm = [mean_OneLRBm;std_OneLRBm];% this makes a two row matrix with mean memory on top and Stddev below. it is returned to the main program.
 
    % extract the left and rigt values of the mean and std promoter levels at the end of
    % the maintenance phase. 
    
    meanLembryoB2 (i) = meanembryoB2 (2*i-1); % extract the left mean values for the promoter
    meanRembryoB2 (i) = meanembryoB2 (2*i); % extract the right mean values for the promoter
    stdLembryoB2 (i)  = stdembryoB2 (2*i-1); % extract the left std for the promoter
    stdRembryoB2 (i)  = stdembryoB2 (2*i); % extract the right std for the promoter
    
    meanstdLembryoB2 = [meanLembryoB2;stdLembryoB2];% this makes a two row matrix with mean on top and Stddev below.
    meanstdRembryoB2 = [meanRembryoB2;stdRembryoB2];% this makes a two row matrix with mean on top and Stddev below.
 
    
end

numR = numel (parameter1); %input the size of the block of parameter 1 values. (will be rows)
numC = numel (parameter2); %input the size of the block of parameter 2 values. (will be columns)

leftBe = reshape (meanLBe, numR,numC); %left side for promoter
rightBe = reshape (meanRBe, numR,numC); %right side for promoter
LRmeanBe = reshape (meanLRBe, numR,numC);%mean of left and right for promoter

leftBm = reshape (meanLBm, numR,numC); %left side for promoter
rightBm = reshape (meanRBm, numR,numC); %right side for promoter
LRmeanBm = reshape (meanLRBm, numR,numC);%mean of left and right for promoter

leftAe = reshape (meanLAe, numR,numC); %left side for promoter
rightAe = reshape (meanRAe, numR,numC); %right side for promoter
LRmeanAe = reshape (meanLRAe, numR,numC);%mean of left and right for promoter

leftAm = reshape (meanLAm, numR,numC); %left side for promoter
rightAm = reshape (meanRAm, numR,numC); %right side for promoter
LRmeanAm = reshape (meanLRAm, numR,numC);%mean of left and right for promoter



lastcol = zeros (numR,1); 
lastrow = zeros (1,numC+1);

leftBe0 = [leftBe,lastcol]; leftBe0 = [leftBe0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
rightBe0 = [rightBe,lastcol]; rightBe0 = [rightBe0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
LRmeanBe0 = [LRmeanBe, LRmeanBe(:,numC)]; LRmeanBe0 = [LRmeanBe0;LRmeanBe0(numR,:)]; % this adds a row of zeros at the last column and row for pcolor.

leftBm0 = [leftBm,lastcol]; leftBm0 = [leftBm0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
rightBm0 = [rightBm,lastcol]; rightBm0 = [rightBm0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
LRmeanBm0 = [LRmeanBm,LRmeanBm(:,numC)]; LRmeanBm0 = [LRmeanBm0;LRmeanBm0(numR,:)]; % this adds a row of zeros at the last column and row for pcolor.

leftAe0 = [leftAe,lastcol]; leftAe0 = [leftAe0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
rightAe0 = [rightAe,lastcol]; rightAe0 = [rightAe0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
LRmeanAe0 = [LRmeanAe,LRmeanAe(:,numC)]; LRmeanAe0 = [LRmeanAe0;LRmeanAe0(numR,:)]; % this adds a row of zeros at the last column and row for pcolor.

leftAm0 = [leftAm,lastcol]; leftAm0 = [leftAm0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
rightAm0 = [rightAm,lastcol]; rightAm0 = [rightAm0;lastrow]; % this adds a row of zeros at the last column and row for pcolor.
LRmeanAm0 = [LRmeanAm,LRmeanAm(:,numC)]; LRmeanAm0 = [LRmeanAm0;LRmeanAm0(numR,:)]; % this adds a row of zeros at the last column and row for pcolor.



% plot fig1: fit of data to controls
%---------------------------------------------------------------------------------------------------------------------------------

lacZ= [1 1 1;0.980392158031464 0.980392158031464 0.980392158031464;0.960784316062927 0.960784316062927 0.960784316062927;0.941176474094391 0.941176474094391 0.941176474094391;0.921568632125854 0.921568632125854 0.921568632125854;0.901960790157318 0.901960790157318 0.901960790157318;0.882352948188782 0.882352948188782 0.882352948188782;0.862745106220245 0.862745106220245 0.862745106220245;0.846666693687439 0.854509830474854 0.876470565795898;0.830588221549988 0.846274495124817 0.890196084976196;0.814509809017181 0.838039219379425 0.903921604156494;0.798431396484375 0.829803943634033 0.917647063732147;0.782352924346924 0.821568608283997 0.9313725233078;0.766274511814117 0.813333332538605 0.945098042488098;0.750196099281311 0.805098056793213 0.958823561668396;0.734117686748505 0.796862781047821 0.972549021244049;0.718039214611053 0.788627445697784 0.986274480819702;0.701960802078247 0.780392169952393 1;0.686274528503418 0.764705896377563 1;0.670588254928589 0.749019622802734 1;0.65490198135376 0.733333349227905 1;0.643137276172638 0.717647075653076 1;0.627451002597809 0.701960802078247 1;0.61176472902298 0.682352960109711 1;0.596078455448151 0.666666686534882 1;0.580392181873322 0.647058844566345 1;0.564705908298492 0.627451002597809 1;0.552941203117371 0.61176472902298 1;0.537254929542542 0.592156887054443 1;0.521568655967712 0.572549045085907 1;0.505882382392883 0.552941203117371 1;0.490196079015732 0.529411792755127 1;0.474509805440903 0.509803950786591 1;0.458823531866074 0.490196079015732 1;0.447058826684952 0.466666668653488 1;0.431372553110123 0.447058826684952 1;0.415686279535294 0.423529416322708 1;0.400000005960464 0.400000005960464 1;0.380392163991928 0.384313732385635 0.980392158031464;0.360784322023392 0.372549027204514 0.960784316062927;0.341176480054855 0.356862753629684 0.941176474094391;0.325490206480026 0.345098048448563 0.917647063732147;0.30588236451149 0.333333343267441 0.898039221763611;0.290196090936661 0.321568638086319 0.878431379795074;0.274509817361832 0.309803932905197 0.858823537826538;0.254901975393295 0.298039227724075 0.839215695858002;0.239215686917305 0.286274522542953 0.819607853889465;0.223529413342476 0.274509817361832 0.796078443527222;0.211764708161354 0.266666680574417 0.776470601558685;0.196078434586525 0.254901975393295 0.756862759590149;0.184313729405403 0.24705882370472 0.74117648601532;0.176470592617989 0.239215686917305 0.725490212440491;0.168627455830574 0.235294118523598 0.713725507259369;0.156862750649452 0.227450981736183 0.69803923368454;0.149019613862038 0.219607844948769 0.682352960109711;0.141176477074623 0.215686276555061 0.666666686534882;0.129411771893501 0.207843139767647 0.65490198135376;0.121568627655506 0.20392157137394 0.639215707778931;0.113725490868092 0.196078434586525 0.623529434204102;0.105882354080677 0.192156866192818 0.607843160629272;0.0980392172932625 0.184313729405403 0.592156887054443;0.0941176488995552 0.180392161011696 0.580392181873322;0.0862745121121407 0.172549024224281 0.564705908298492;0.0784313753247261 0.168627455830574 0.549019634723663];

figure 

colormap (flipud (parula));

 title1 = {'Plot option 2: parameter space: silencing';note;'';'promoter memory'};
 title2 = {'Plot option 2: parameter space: silencing';note;'';'promoter initiation: all data'};
 

 
 subplot1 = subplot(1,2,1);%late embryo promoter compared to control
 pcolor (LRmeanBm0); caxis ([0,0.3]);shading interp;
 title (title1);
 set(gca,'XTick',1:1:num2+1)
 set(gca,'XTickLabel',(parameter2_list));
 set(gca,'YTick',1:1:num1+1)
 set(gca,'YTickLabel',(parameter1_list));
 xlabel(parameter2_label);
 ylabel(parameter1_label);

 
 subplot2 = subplot(1,2,2);%late embryo PRE/TRE compared to control
 pcolor (LRmeanAm0); caxis ([0,0.3]);shading interp;
 title ('PRE/TRE memory');
 set(gca,'XTick',1:1:num2+1)
 set(gca,'XTickLabel',(parameter2_list));
 set(gca,'YTick',1:1:num1+1)
 set(gca,'YTickLabel',(parameter1_list));
 xlabel(parameter2_label);
 ylabel(parameter1_label);


figure 
colormap (lacZ)

lastcol2 = zeros (n,1);
embryoB10 = [embryoB1,lastcol2];
embryoB20 = [embryoB2,lastcol2];
embryoA10 = [embryoA1,lastcol2];
embryoA20 = [embryoA2,lastcol2];

subplot3 = subplot(2,2,1); 
 pcolor (embryoB10);caxis ([0,1]); shading flat;
 title (title2);
 set(gca,'XTick',1:num1*2:(num2*num1*2)+num1*2);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('repetitions');
 
 subplot4 = subplot(2,2,3);
 pcolor (embryoB20);caxis ([0,1]); shading flat;
 title ('promoter maintenance: all data');
 set(gca,'XTick',1:num1*2:(num2*num1*2)+num1*2);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('repetitions');
 
 subplot5 = subplot(2,2,2);
 pcolor (embryoA10); caxis ([0,0.7]);shading flat;
 title ('PRE/TRE initiation: all data');
 set(gca,'XTick',1:num1*2:(num2*num1*2)+num1*2);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('repetitions');
 
 subplot6 = subplot(2,2,4);
 pcolor (embryoA20); caxis ([0,0.7]);shading flat;
 title ('PRE/TRE maintenance: all data');
 set(gca,'XTick',1:num1*2:(num2*num1*2)+num1*2);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('repetitions');
 
 
figure

title3 = {'Plot option 2: parameter space: silencing';note;'';'Promoter memory: means (1= good, 0.4 = bad)'};

subplot7 = subplot(3,1,1);
hold on
ax1 = subplot(3,1,1);
bar(singlevalBm (1,:))
errorbar(singlevalBm (1,:),singlevalBm (2,:),'.k');
ylim(ax1,[0.4 1.2])
title (title3);
 set(gca,'XTick',1:num1:(num2*num1)+num1);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('silencing memory');
 
subplot8 = subplot(3,1,2);
hold on
ax2 = subplot(3,1,2);
bar(meanstdLembryoB2 (1,:))
errorbar(meanstdLembryoB2 (1,:),meanstdLembryoB2 (2,:),'.k');
ylim(ax2,[0 1.6])
title ('Promoter levels at end of maintenance phase: Anterior: should be silenced');
 set(gca,'XTick',1:num1:(num2*num1)+num1);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('N_B/s');

subplot9 = subplot(3,1,3);
hold on
ax3 = subplot(3,1,3);
bar(meanstdRembryoB2 (1,:))
errorbar(meanstdRembryoB2 (1,:),meanstdRembryoB2 (2,:),'.k');
ylim(ax3,[0 1.6])
title ('Promoter levels at end of maintenance phase: Posterior: should be active');
 set(gca,'XTick',1:num1:(num2*num1)+num1);
 set(gca,'XTickLabel',(parameter2_list));
 xlabel(Xlabel2);
 ylabel('N_B/s');
 

end

