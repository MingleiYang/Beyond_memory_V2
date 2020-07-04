function [Bend, Aend, AMend,meanBend,stdBend, meanAend,stdAend,meanAMend,stdAMend,...
     memory_score_nw, mean_memory_nw, std_memory_nw]...
     = plot_option_1_3(note, plot_option,BIG_meanA,BIG_meanM,BIG_meanB,n,numstepeq,numstep1,numstep2,numstep3,numstep)

% plot_option_1_3  Memory of silencing and activation: embryo time course
% up to 7.5 hours of development. Used in Figures 2,3 S2,S3 of Reinig et
% al., 2020
% Leonie Ringrose, 17.05.19, modified to include box and whisker plots,
% 22.06.2020. 

% extract data from output BIG_meanA, BIG_meanM and BIG_meanB for plots. 
%------------------------------------------------------------------------------------------------------------------
    BIG_meanAM = BIG_meanA - BIG_meanM; % This converts the scale for the PRE/TRE to [-1 : +1]

    numstep1 = numstepeq + numstep1;

    Brep1 = BIG_meanB (1:numstep1, :);
    Brep2 = BIG_meanB (numstep1+1:numstep-numstep3, :);
    Brep3 = BIG_meanB (numstep1+numstep2+1:numstep, :);
    Bend  = BIG_meanB (numstep, :); 

    AMrep1 = BIG_meanAM (1:numstep1, :);
    AMrep2 = BIG_meanAM (numstep1+1:numstep-numstep3, :);
    AMrep3 = BIG_meanAM (numstep1+numstep2+1:numstep, :);
    AMend  = BIG_meanAM (numstep, :); 
    Aend  = BIG_meanA (numstep, :); 

    x = numstep; 


    for i = 1:n-1
   
  
        Brep1 = [Brep1;(BIG_meanB((i*x)+1:(numstep1+(i*x)),:))];
        Brep2 = [Brep2; (BIG_meanB (numstep1+(i*x)+1:numstep-numstep3+(i*x), :))];
        Brep3 = [Brep3; (BIG_meanB (numstep1+numstep2+(i*x)+1:numstep +(i*x), :))];
        Bend = [Bend;(BIG_meanB ((i+1)*x,:))]; 
    
 
        AMrep1 = [AMrep1;(BIG_meanAM((i*x)+1:(numstep1+(i*x)),:))];
        AMrep2 = [AMrep2; (BIG_meanAM (numstep1+(i*x)+1:numstep-numstep3+(i*x), :))];
        AMrep3 = [AMrep3; (BIG_meanAM (numstep1+numstep2+(i*x)+1:numstep +(i*x), :))];
        AMend = [AMend;(BIG_meanAM ((i+1)*x,:))];
        Aend = [Aend;(BIG_meanA ((i+1)*x,:))];
    
    end
    
    %calculate means and std of last 10 mins for plot:
    
    meanBend = mean (Bend);
    stdBend = std (Bend);
    
    meanAMend = mean (AMend);
    stdAMend = std (AMend);
    
    meanAend = mean (Aend);
    stdAend = std (Aend);
    
    %make a row of zeros of size n
    lastrow = zeros (1,2*n);
    % add it to each "embryo" for pcolor
   
    % Promoter, no coupling. 
    %---------------------------------------------------------------------------------------------------------------------
    
    Brep1noL = reshape (Brep1 (:,1), numstep1,n);
    Brep2noL = reshape (Brep2 (:,1), numstep2,n);
    Brep3noL = reshape (Brep3 (:,1), numstep3,n);

    Brep1noR = reshape (Brep1 (:,2), numstep1,n);
    Brep2noR = reshape (Brep2 (:,2), numstep2,n); 
    Brep3noR = reshape (Brep3 (:,2), numstep3,n); 
    
    
    Brep1noLR = [flipud([Brep1noL,Brep1noR]);lastrow]; % Rep1 combine left and right to a single "embryo", flip vertically, add zeros for pcolor. 
    lastcol1 = zeros (numstep1+1,1);                    % in the matrix, time windows run vertically from bottom to top; Last row is zeros.
                                                        %in the pcolor plot they will appear vertically from top to bottom, last row no plotted. repetitions across. 
    Brep1noLR = [Brep1noLR,lastcol1];                   % Add column of zeros for pcolor. 
    
    Brep2noLR = [flipud([Brep2noL,Brep2noR]);lastrow]; %Rep 2 combine left and right. 
    lastcol2 = zeros (numstep2+1,1);            
    Brep2noLR = [Brep2noLR,lastcol2];
    
    Brep3noLR = [flipud([Brep3noL,Brep3noR]);lastrow]; %Rep 3 combine left and right. 
    lastcol3 = zeros (numstep3+1,1);            
    Brep3noLR = [Brep3noLR,lastcol3];
    
    % Make single time course out of all three blocks (this is not plotted in the output but it may be useful). 
    %---------------------------------------------------------------------------------------------------------------------
    lastrow0 = zeros (1,2*n+1);
    Brep1_3noLR = [(Brep3noLR (1:numstep3,:));(Brep2noLR (1:numstep2,:));(Brep1noLR (1:numstep1,:))];
    Brep1_3noLR = [Brep1_3noLR;lastrow0];
    
    % Promoter, with coupling. 
    %---------------------------------------------------------------------------------------------------------------------
    Brep1withL = reshape (Brep1 (:,3), numstep1,n);
    Brep2withL = reshape (Brep2 (:,3), numstep2,n);
    Brep3withL = reshape (Brep3 (:,3), numstep3,n);
    
    Brep1withR = reshape (Brep1 (:,4), numstep1,n);
    Brep2withR = reshape (Brep2 (:,4), numstep2,n);
    Brep3withR = reshape (Brep3 (:,4), numstep3,n);
    
    Brep1withLR = [flipud([Brep1withL,Brep1withR]);lastrow]; % Promoter with coupling Rep1 combine left and right. 
    Brep1withLR = [Brep1withLR,lastcol1];
    
    Brep2withLR = [flipud([Brep2withL,Brep2withR]);lastrow]; % Promoter with coupling Rep2 combine left and right. 
    Brep2withLR = [Brep2withLR,lastcol2];
    
    Brep3withLR = [flipud([Brep3withL,Brep3withR]);lastrow]; % Promoter with coupling Rep3 combine left and right. 
    Brep3withLR = [Brep3withLR,lastcol3];
    
    % Make single time course out of all three blocks (this is not plotted in the output but it may be useful).
    %---------------------------------------------------------------------------------------------------------------------
   
    Brep1_3withLR = [(Brep3withLR (1:numstep3,:));(Brep2withLR (1:numstep2,:));(Brep1withLR (1:numstep1,:))];
    Brep1_3withLR = [Brep1_3withLR;lastrow0];
    
    %PRE/TRE, no coupling. 
    %---------------------------------------------------------------------------------------------------------------------
    
    AMrep1noL = reshape (AMrep1 (:,1), numstep1,n);
    AMrep2noL = reshape (AMrep2 (:,1), numstep2,n);
    AMrep3noL = reshape (AMrep3 (:,1), numstep3,n);
    
    AMrep1noR = reshape (AMrep1 (:,2), numstep1,n);
    AMrep2noR = reshape (AMrep2 (:,2), numstep2,n);
    AMrep3noR = reshape (AMrep3 (:,2), numstep3,n);
    
    AMrep1noLR = [flipud([AMrep1noL,AMrep1noR]);lastrow]; % PRE/TRE without coupling Rep1 combine left and right. 
    AMrep1noLR = [AMrep1noLR,lastcol1];
    
    AMrep2noLR = [flipud([AMrep2noL,AMrep2noR]);lastrow]; % PRE/TRE without coupling Rep2 combine left and right. 
    AMrep2noLR = [AMrep2noLR,lastcol2];
    
    AMrep3noLR = [flipud([AMrep3noL,AMrep3noR]);lastrow]; % PRE/TRE without coupling Rep3 combine left and right. 
    AMrep3noLR = [AMrep3noLR,lastcol3];
    
    % Make single time course out of all three blocks (this is not plotted in the output but it may be useful).
    %---------------------------------------------------------------------------------------------------------------------
  
    AMrep1_3noLR = [(AMrep3noLR (1:numstep3,:));(AMrep2noLR (1:numstep2,:));(AMrep1noLR (1:numstep1,:))];
    AMrep1_3noLR = [AMrep1_3noLR;lastrow0];
    
    
    %PRE/TRE, with coupling. 
    %---------------------------------------------------------------------------------------------------------------------
    
    AMrep1withL = reshape (AMrep1 (:,3), numstep1,n);
    AMrep2withL = reshape (AMrep2 (:,3), numstep2,n);
    AMrep3withL = reshape (AMrep3 (:,3), numstep3,n);
    
    AMrep1withR = reshape (AMrep1 (:,4), numstep1,n);
    AMrep2withR = reshape (AMrep2 (:,4), numstep2,n);
    AMrep3withR = reshape (AMrep3 (:,4), numstep3,n);
     
    AMrep1withLR = [flipud([AMrep1withL,AMrep1withR]);lastrow]; % PRE/TRE without coupling Rep1 combine left and right. 
    AMrep1withLR = [AMrep1withLR,lastcol1];
    
    AMrep2withLR = [flipud([AMrep2withL,AMrep2withR]);lastrow]; % PRE/TRE without coupling Rep2 combine left and right. 
    AMrep2withLR = [AMrep2withLR,lastcol2];
    
    AMrep3withLR = [flipud([AMrep3withL,AMrep3withR]);lastrow]; % PRE/TRE without coupling Rep3 combine left and right. 
    AMrep3withLR = [AMrep3withLR,lastcol3];
    
    % Calculate memory scores (applies to promoter only, in the anterior)
    %--------------------------------------------------------------------
    
    diff_n = Brep3noLR(1,1:n) - Brep2noLR (1,1:n); 
    memory_score_n = 1 - diff_n.^2;
    mean_memory_n = mean (memory_score_n); 
    std_memory_n = std (memory_score_n);
    
    diff_w = Brep3withLR(1,1:n) - Brep2noLR (1,1:n); 
    memory_score_w = 1 - diff_w.^2;
    mean_memory_w = mean (memory_score_w); 
    std_memory_w = std (memory_score_w);
    
    memory_score_nw = [memory_score_n; memory_score_w];
    mean_memory_nw = [mean_memory_n; mean_memory_w];
    std_memory_nw = [std_memory_n; std_memory_w];
    
    
    % Make single time course out of all three blocks (this is not plotted in the output but it may be useful).
    %---------------------------------------------------------------------------------------------------------------------
      
    AMrep1_3withLR = [(AMrep3withLR (1:numstep3,:));(AMrep2withLR (1:numstep2,:));(AMrep1withLR (1:numstep1,:))];
    AMrep1_3withLR = [AMrep1_3withLR;lastrow0];
    
    
    lacZ= [1 1 1;0.980392158031464 0.980392158031464 0.980392158031464;0.960784316062927 0.960784316062927 0.960784316062927;0.941176474094391 0.941176474094391 0.941176474094391;0.921568632125854 0.921568632125854 0.921568632125854;0.901960790157318 0.901960790157318 0.901960790157318;0.882352948188782 0.882352948188782 0.882352948188782;0.862745106220245 0.862745106220245 0.862745106220245;0.846666693687439 0.854509830474854 0.876470565795898;0.830588221549988 0.846274495124817 0.890196084976196;0.814509809017181 0.838039219379425 0.903921604156494;0.798431396484375 0.829803943634033 0.917647063732147;0.782352924346924 0.821568608283997 0.9313725233078;0.766274511814117 0.813333332538605 0.945098042488098;0.750196099281311 0.805098056793213 0.958823561668396;0.734117686748505 0.796862781047821 0.972549021244049;0.718039214611053 0.788627445697784 0.986274480819702;0.701960802078247 0.780392169952393 1;0.686274528503418 0.764705896377563 1;0.670588254928589 0.749019622802734 1;0.65490198135376 0.733333349227905 1;0.643137276172638 0.717647075653076 1;0.627451002597809 0.701960802078247 1;0.61176472902298 0.682352960109711 1;0.596078455448151 0.666666686534882 1;0.580392181873322 0.647058844566345 1;0.564705908298492 0.627451002597809 1;0.552941203117371 0.61176472902298 1;0.537254929542542 0.592156887054443 1;0.521568655967712 0.572549045085907 1;0.505882382392883 0.552941203117371 1;0.490196079015732 0.529411792755127 1;0.474509805440903 0.509803950786591 1;0.458823531866074 0.490196079015732 1;0.447058826684952 0.466666668653488 1;0.431372553110123 0.447058826684952 1;0.415686279535294 0.423529416322708 1;0.400000005960464 0.400000005960464 1;0.380392163991928 0.384313732385635 0.980392158031464;0.360784322023392 0.372549027204514 0.960784316062927;0.341176480054855 0.356862753629684 0.941176474094391;0.325490206480026 0.345098048448563 0.917647063732147;0.30588236451149 0.333333343267441 0.898039221763611;0.290196090936661 0.321568638086319 0.878431379795074;0.274509817361832 0.309803932905197 0.858823537826538;0.254901975393295 0.298039227724075 0.839215695858002;0.239215686917305 0.286274522542953 0.819607853889465;0.223529413342476 0.274509817361832 0.796078443527222;0.211764708161354 0.266666680574417 0.776470601558685;0.196078434586525 0.254901975393295 0.756862759590149;0.184313729405403 0.24705882370472 0.74117648601532;0.176470592617989 0.239215686917305 0.725490212440491;0.168627455830574 0.235294118523598 0.713725507259369;0.156862750649452 0.227450981736183 0.69803923368454;0.149019613862038 0.219607844948769 0.682352960109711;0.141176477074623 0.215686276555061 0.666666686534882;0.129411771893501 0.207843139767647 0.65490198135376;0.121568627655506 0.20392157137394 0.639215707778931;0.113725490868092 0.196078434586525 0.623529434204102;0.105882354080677 0.192156866192818 0.607843160629272;0.0980392172932625 0.184313729405403 0.592156887054443;0.0941176488995552 0.180392161011696 0.580392181873322;0.0862745121121407 0.172549024224281 0.564705908298492;0.0784313753247261 0.168627455830574 0.549019634723663];

    figure
    colormap (lacZ)

    if plot_option == 1
        title1 = {'Plot option 1: embryo time course';note;'';'Promoter : no coupling';'Cycles 1-13'};
    elseif plot_option == 3
        title1 = {'Plot option 3: embryo time course, heat shock';note;'';'Promoter : no coupling';'Cycles 1-13'};
    end  
    
    title2 = {'Promoter : with coupling'; 'Cycles 1-13'};
    title3 = {'PRE/TRE : no coupling'; 'Cycles 1-13'};
    title4 = {'PRE/TRE : with coupling'; 'Cycles 1-13'};
    
% Plots the time course data in separate blocks as in the paper.
% row 1: cycles 1-13.
%---------------------------------------------------------------------------------------------------------------------

  subplot1 = subplot (3,4,1);
  pcolor (Brep1noLR); caxis ([0,1]);shading flat;   % 1) promoter no coupling
  title (title1);
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  ylabel('time (0 - 2h10)');
  
  subplot2 = subplot (3,4,2);
  pcolor (Brep1withLR); caxis ([0,1]);shading flat; % 2) promoter with coupling
  title (title2);
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  
  subplot3 = subplot (3,4,3);
  pcolor (AMrep1noLR); caxis ([-1,1]);shading flat; % 3) PRE/TRE no coupling 
  title (title3);
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  
  subplot4 = subplot (3,4,4);
  pcolor (AMrep1withLR); caxis ([-1,1]);shading flat; % 4) PRE/TRE with coupling 
  title (title4);
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  
% row 2: cycle 14.
%---------------------------------------------------------------------------------------------------------------------

  subplot5 = subplot (3,4,5);
  pcolor (Brep2noLR); caxis ([0,1]);shading flat; % 1) promoter no coupling
  title ('cycle 14');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  ylabel('time (2h10-4h50)');
  
  subplot6 = subplot (3,4,6);
  pcolor (Brep2withLR); caxis ([0,1]);shading flat;% 2) promoter with coupling
  title ('cycle 14');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  
  subplot7 = subplot (3,4,7);
  pcolor (AMrep2noLR); caxis ([-1,1]);shading flat;% 3) PRE/TRE no coupling
  title ('cycle 14');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  
  subplot8 = subplot (3,4,8);
  pcolor (AMrep2withLR); caxis ([-1,1]);shading flat;% 4) PRE/TRE with coupling 
  title ('cycle 14');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  

% row 3: cycle 15.
%---------------------------------------------------------------------------------------------------------------------

  subplot9 = subplot (3,4,9);
  pcolor (Brep3noLR); caxis ([0,1]);shading flat;% 1) promoter no coupling
  title ('cycle 15');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  xlabel('anterior    |     posterior');
  ylabel('time (4h50-7h30)');
  
  subplot10 = subplot (3,4,10);
  pcolor (Brep3withLR); caxis ([0,1]);shading flat;% 2) promoter with coupling
  title ('cycle 15');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  xlabel('anterior    |     posterior');
  
  subplot11 = subplot (3,4,11);
  pcolor (AMrep3noLR); caxis ([-1,1]);shading flat;% 3) PRE/TRE no coupling
  title ('cycle 15');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  xlabel('anterior    |     posterior');
  
  subplot12 = subplot (3,4,12);
  pcolor (AMrep3withLR); caxis ([-1,1]);shading flat;% 4) PRE/TRE with coupling 
  title ('cycle 15');
  set(gca,'YTickLabel',{' '});
  set(gca,'XTickLabel',{' '});
  xlabel('anterior    |     posterior');
  

figure

% Plot bar plots and boxplots of data for last 10 mins
if plot_option == 1
    title5 = {'Plot option 1: embryo time course';note;'';'Promoter B state: mean of last 10 mins'};

elseif plot_option == 3
    title5 = {'Plot option 3: embryo time course';note;'';'Promoter B state: mean of last 10 mins'};
end

    subplot13 = subplot(2,2,1);
    hold on
    bar(meanBend); 
    errorbar(meanBend,stdBend,'.k'); 
    ylim ([0 1.6]);
    title (title5);
    set(gca,'XTickLabel',{' ','A-no coupling',' ','P-no coupling',' ','A-with coupling',' ','P-with coupling '});
    ylabel('Mean promoter (N_B/s)');
    
    subplot14 = subplot(2,2,2);
    hold on
    bar(meanAend); % for A and M states on scale of -1: +1, plot meanAMend, stdAMend.
    errorbar(meanAend,stdAend,'.k'); 
    ylim ([0 1.6]);
    title ('PRE/TRE A state: mean of last 10 mins');
    set(gca,'XTickLabel',{' ','A-no coupling',' ','P-no coupling',' ','A-with coupling',' ','P-with coupling '});
    ylabel('Mean PRE/TRE (N_A/S)');

    subplot15 = subplot(2,2,3);
    hold on
    boxplot(Bend,'Symbol','ko','Colors','k','OutlierSize',3);  
    %ylim ([0 1.6]);
    title ('Promoter: boxplot of last 10 mins');
    set(gca,'XTickLabel',{'A-no coupling','P-no coupling','A-with coupling','P-with coupling '});
    ylabel('Mean promoter (N_B/s)');
    
    subplot16 = subplot(2,2,4);
    hold on
    boxplot(AMend,'Symbol','ko','Colors','k','OutlierSize',3); % for A and M states on scale of -1: +1, plot meanAMend, stdAMend.
    %ylim ([0 1.6]);
    title ('PRE/TRE: boxplot of last 10 mins');
    set(gca,'XTickLabel',{'A-no coupling','P-no coupling','A-with coupling','P-with coupling '});
    ylabel('Mean PRE/TRE (N_A-N_M/S)');       
  
    if plot_option == 1
        
        figure
        %plot bar and boxplots for memory scores
        %------------------------------------------
    
        subplot17 = subplot(2,1,1);
        hold on
        bar(mean_memory_nw); 
        errorbar(mean_memory_nw,std_memory_nw,'.k'); 
        %ylim ([0 1.6]);
        title ('Memory score: Anterior');
        set(gca,'XTickLabel',{' ','no coupling',' ','with coupling'});
        ylabel('Memory score');
    
        subplot18 = subplot(2,1,2);
        hold on
        boxplot(memory_score_nw','Symbol','ko','Colors','k','OutlierSize',3);  
        %ylim ([0 1.6]);
        title ('Memory score: Anterior');
        set(gca,'XTickLabel',{'no coupling','with coupling'});
        ylabel('Memory score');
          
    end
    
    
end

