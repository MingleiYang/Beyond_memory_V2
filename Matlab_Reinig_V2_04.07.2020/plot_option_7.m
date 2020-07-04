function [figure1] = plot_option_7(t,k,n,sum_timeA,sum_timeB,sum_timeM,sum_timeF,...
    ALL_stateP_no,ALL_stateP_with,ALL_stateE_no,ALL_stateE_with,note,heatshock)

% plot_option_7   Plots mean time course for PRE/TRE and promoter, and gives state analysis for last 10 mins of simulation. 
% Leonie Ringrose, 17.05.19


% Convert data to means

mean_timeA = sum_timeA/k;
mean_timeM = sum_timeM/k;
mean_timeB = sum_timeB/k;
mean_timeF = sum_timeF/k;

columnO = 1;%gives the timecourse for the first p1 value without coupling
columnC = 3;%gives the timecourse for the first p1 value with coupling.


% Convert time from iterations to hours, extract values for X axis.

timelist = '0';

for Str = 300:300:t 
 
    timen = mat2str(Str/300);
    timelist = char(timelist,timen);
      
end

timelist2 = '0';

for Str = 300*24:300*24:t 
 
    timen2 = mat2str(Str/(24*300));
    timelist2 = char(timelist2,timen2);
 
end

    if t >= 11520
        timelegend = timelist2;
        Xtick = 1:300*24:t;
        time = ('time (days)');
    else
        timelegend = timelist;
        Xtick = 1:300:t;
        time = ('time (hours)');
    end


    if heatshock == 1
        title1 = {'Plot option 7, embryo time course: early heat shock';note;'';'PRE/TRE: no coupling'};     
    elseif heatshock == 2
        title1 = {'Plot option 7, embryo time course: late heat shock';note;'';'PRE/TRE: no coupling'};
    else
        title1 = {'Plot option 7, embryo time course: no heat shock';note;'';'PRE/TRE: no coupling'};
    end
    
figure1 = figure

subplot1 = subplot(2,2,1); 
    plot(mean_timeA (:,columnO)); 
    hold all; 
    plot (mean_timeM (:,columnO)); 
    set(gca,'XTick',(Xtick));
    set(gca,'XTickLabel',(timelegend));
    ylim ([0 1]);
    xlim ([0, t]);
    title (title1);
    legend('A', 'M');
    xlabel(time);
    ylabel('mean A or M');

subplot2 = subplot(2,2,2); 
    plot(mean_timeA (:,columnC)); 
    hold all; 
    plot (mean_timeM (:,columnC)); 
    set(gca,'XTick',(Xtick));
    set(gca,'XTickLabel',(timelegend));
    ylim ([0 1]);
    xlim ([0, t]);
    title ('PRE/TRE: with coupling');
    legend('A', 'M');
    xlabel(time);
    ylabel('mean A or M');
    
subplot3 = subplot(2,2,3); 
    plot(mean_timeB (:,columnO)); 
    hold all; 
    plot (mean_timeF (:,columnO)); 
    set(gca,'XTick',(Xtick));
    set(gca,'XTickLabel',(timelegend));
    ylim ([0 1.2]);
    xlim ([0, t]);
    title ('Promoter: no coupling');
    legend('B', 'F');
    xlabel(time);
    ylabel('mean B or F');

subplot4 = subplot(2,2,4); 
    plot(mean_timeB (:,columnC)); 
    hold all; 
    plot (mean_timeF (:,columnC)); 
    set(gca,'XTick',(Xtick));
    set(gca,'XTickLabel',(timelegend));
    ylim ([0 1.2]);
    xlim ([0, t]);
    title ('Promoter: with coupling');
    legend('B', 'F');
    xlabel(time);
    ylabel('mean B or F');







A = size (ALL_stateP_no);
last = A (1);

%calculate the % of each state represented, over the last time window

statePN = ALL_stateP_no(last,:);
statePW = ALL_stateP_with(last,:);
stateEN = ALL_stateE_no(last,:);
stateEW = ALL_stateE_with(last,:);


countsPN = zeros (n,3); %set up output matrix.
countsPW = zeros (n,3); %set up output matrix.
countsEN = zeros (n,3); %set up output matrix.
countsEW = zeros (n,3); %set up output matrix.


for h = 1:n; % count the states for each iteration
    if statePN (h) == 1  % A
    countsPN (h,1)= 1;
    elseif statePN (h)== -1;
    countsPN (h,3)= 1; % M
    else countsPN (h,2)=1; % U
    end 
    
    if statePW (h) == 1  % A
    countsPW (h,1)= 1;
    elseif statePW (h)== -1;% M
    countsPW (h,3)= 1; 
    else countsPW (h,2)=1; % U
    end 
    
    if stateEN (h) == 1  % B
    countsEN (h,1)= 1;
    elseif stateEN (h)== -1;% F
    countsEN (h,3)= 1; 
    else countsEN (h,2)=1; % N
    end 
    
    if stateEW (h) == 1  % B
    countsEW (h,1)= 1;
    elseif stateEW (h)== -1;% F
    countsEW (h,3)= 1; 
    else countsEW (h,2)=1; % N
    end 
end 

% sum up the counts in each category. 

sumPN = sum (countsPN)+0.0001; % gives three values, number of A, M and U states. 
sumPW = sum (countsPW)+0.0001; % dirty hack to avoid zeros for pie chart!
sumEN = sum (countsEN)+0.0001;
sumEW = sum (countsEW)+0.0001;

    if heatshock == 1
        title2 = {'Plot option 7, embryo time course: early heat shock';note;'';'PRE/TRE - no coupling: final states'};     
    elseif heatshock == 2
        title2 = {'Plot option 7, embryo time course: late heat shock';note;'';'PRE/TRE - no coupling: final states'};
    else
        title2 = {'Plot option 7, embryo time course: no heat shock';note;'';'PRE/TRE - no coupling: final states'};
    end


figure; 
colormap (autumn);
subplot5 = subplot(2,2,1);
pie (sumPN);
title (title2);
legend('A', 'U', 'M');

subplot6 = subplot(2,2,2);
pie (sumPW);
title ('PRE/TRE - with coupling: final states');
legend('A', 'U', 'M');

subplot7 = subplot(2,2,3);
pie (sumEN);
title ('Promoter - no coupling: final states');
legend('B', 'N', 'F');

subplot8 = subplot(2,2,4);
pie (sumEW);
title ('Promoter - with coupling: final states');
legend('B', 'N', 'F');

end

