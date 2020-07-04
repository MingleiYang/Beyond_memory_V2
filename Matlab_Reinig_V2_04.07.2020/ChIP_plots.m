% Leonie Ringrose, 25.06.2020. 
% This is a standalone script accompanying the manuscript Reinig et al, Nature Communications 2020
% The code contains raw data and plotting function for the final version of Figure 4E (Bar plots
% with standard deviation, overlaid with individual data points). Revised
% upon acceptance of the paper to comply with editorial request. 
% The data contain chip enrichments as a percentage of H3 ChIP for 5 time
% points and 6 antibodies, performed in Drosophila embryos, each row
% contain 4 replicates for each time point. 0-2h, 2-4h,4-6h,'6-8h','8-10h'


PC=1; K27me3 =2; K4me1= 3; K4me3 = 4; K36me2 = 5; K36me3 = 6; % Names are included in the table for easier identification of data. 

chip_data = [PC    0.007	0.010	0.007	0.009	0.008	0.018	0.016	0.008	0.031	0.063	0.058	0.040	0.112	0.109	0.119	0.078	0.121	0.050	0.060	0.119
K27me3   0.028	0.024	0.016	0.034	0.020	0.030	0.030	0.021	0.032	0.053	0.041	0.038	0.042	0.031	0.031	0.048	0.136	0.175	0.198	0.121
K4me1	 0.154	0.319	0.232	0.226	0.238	0.243	0.232	0.255	0.185	0.546	0.445	0.224	0.694	0.281	0.429	0.638	0.593	0.375	0.433	0.518
K4me3	 0.057	0.065	0.061	0.128	0.038	0.056	0.067	0.037	0.042	0.065	0.042	0.049	0.184	0.158	0.120	0.150	0.113	0.223	0.254	0.081
K36me2   0.044	0.075	0.042	0.076	0.041	0.061	0.063	0.042	0.079	0.025	0.020	0.092	0.078	0.110	0.152	0.117	0.114	0.276	0.306	0.094
K36me3   0.027	0.070	0.050	0.039	0.041	0.075	0.064	0.045	0.008	0.008	0.009	0.009	0.037	0.115	0.124	0.097	0.144	0.248	0.276	0.128];


means = zeros (6,5);
stds = zeros (6,5);

lo = [2 6 10 14 18];
hi = [5 9 13 17 21];


for i = 1:6
    for j = 1:5
    means (i,j)= mean(chip_data (i,lo(j):hi(j)));
    stds (i,j)= std(chip_data (i,lo(j):hi(j)));
    end 
end

titles = {'PC','K27me3','K4me1','K4me3','K36me2','K36me3'};

figure

x = [0.9 0.95 1.05 1.1 1.9 1.95 2.05 2.1 2.9 2.95 3.05 3.1 3.9 3.95 4.05 4.1 4.9 4.95 5.05 5.1];


for i = 1:6
    subplot_i = subplot(6,1,(i));
    hold on
    bar(means(i,:), 'FaceColor',[.5 .5 .5]); 
    errorbar(means(i,:),stds(i,:),'.k'); 
    scatter (x,chip_data(i,2:21),25, 'filled','k');
    %ylim ([0 1.6]);
    title (titles (i));
    set(gca,'XTickLabel',{'','0-2h ','','2-4h','','4-6h ','','6-8h','','8-10h'});
    ylabel('% H3 ChIP');
end