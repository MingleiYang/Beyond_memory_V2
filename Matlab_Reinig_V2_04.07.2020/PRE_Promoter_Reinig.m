
% PRE_Promoter_Reinig 

% Copyright Leonie Ringrose, 04.07.2020

   %Licensed under the Apache License, Version 2.0 (the "License");
   %you may not use this file except in compliance with the License.
   %You may obtain a copy of the License at

   %http://www.apache.org/licenses/LICENSE-2.0

   %Unless required by applicable law or agreed to in writing, software
   %distributed under the License is distributed on an "AS IS" BASIS,
   %WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   %See the License for the specific language governing permissions and
   %limitations under the License.
 
% This script and the accompanying functions run all the simulations in Reinig et al.  
% "A theoretical model of Polycomb/Trithorax action unites stable epigenetic memory and dynamic regulation"
% User defined plot_options 1-7 correspond to the figures as listed below. 

% Plot options
%---------------------------------------------------------------------------------------

% 1: embryo: memory of silencing. Used for plots in Fig. 2 and S2 of the paper.
% 2: embryo: memory of silencing, parameter space analysis. Used for plots in Fig S1 and S2.
% 3: embryo: memory of activation, heat shock experiment. Used for plots in Fig. 3 and S3.
% 4: embryo: memory of activation, parameter space analysis. Used for plots in Fig S1 and S3.
% 5: eye disc: used for plots in Fig. 5 and S8. 
% 6: eye disc: parameter space analysis. Used for plots in Fig. S1 and S8 of the paper. 
    % N.B. this is a long simulation: n= 100 requires 1h approx!
% 7: general: time course for whole run and state analysis for the last 10 min window. 
    % Used for histone level predictions(Fig. 4a)
 

clear all;

% User-defined inputs. To run in default mode, just select a plot_option (1-7). 
% Use defaults for all other inputs. 
% This will run different simulations corresponding to the figures in the paper. 
%-------------------------------------------------------------------------------------------------------------------------------
note = '20.06.20, Test'; 	% Use this space for any notes. 
plot_option = 3;            % Choose 1-7. See above for list. 
model = 1;                  % Define model: Choose 1 - 14. Default = 1. Model 2 gives similar results.                  
n = 50;                     % Number of total repeats of simulation. Default = 50. 100 for parameter scans (see user manual).

Coupling_Strength = 1;      % User defined coupling factor. Default = 1. 
                            % NB the coupling factors given here will reproduce the figures in the User Manual. The revised   
                            % version of the paper (June 2020) uses a coupling strength of 0.5 for Figures 2 and 3 only.  
PRE_init = ('U');           % Define the initial state of the PRE/TRE (A, U or M). Default = U
Enhancer_init = ('F');      % Define the initial state of the promoter (B or F). Default = F

heatshock = 1;     % Plot option 3, 4 and 7. Choose 0 (none) 1 (early) or 2 (late). Ignored by other plot options. 
PRE_ID = ('b');    % Plot option 5 and 7. PRE/TRE identity (e = eya or b = bxd). All other plot options use bxd parameters. 

PcG = 1;           % Plot option 1, 3, 5, 7. Choose 1 (wild type) or 0 (mutant). Loss of PcG at onset of cycle 15  
TrxG = 1;          % Plot option 1, 3, 5, 7. Choose 1 (wild type) or 0 (mutant). Loss of TrxG at onset of cycle 15 

%End of user-defined section. Make changes beyond this point at your peril!
%--------------------------------------------------------------------------------------------------------------------------------

% Define system size. 
%--------------------------------------------------------------------------------------------------------

m = 40;             %number of nucleosomes in the PRE/TRE (default = 40)
tf = 10;            %number of transcription factor binding sites in the promoter (default = 10)

% Define p5 for plot_option 5 according to PRE/TRE identity . 
%--------------------------------------------------------------------------------------------------------
if PRE_ID  == ('e')
    p5in = 0.17;        % p5in = 0.17 for eya PRE/TRE
elseif PRE_ID  == ('b')
    p5in = 0.04;    % p5in = 0.04 for bxd PRE/TRE
else error 'please enter PRE_ID as either e or b'
end

% Coupling regimes are linked automatically to model and plot option
%--------------------------------------------------------------------

if model == 1 || model >= 3 && model <=8
    CF = 1 * Coupling_Strength;

elseif model == 2 || model >= 9 && model <=14
    CF = 0.5 * Coupling_Strength;

else error 'please enter model number between 1 and 14'
end
    

if plot_option <= 4 || plot_option == 7 % embryo
    C1 = 0 * CF;             % coupling strength in trun1 (early cleavage cycles). 
    C2 = 8 * CF;             % coupling strength in trun2 (initiation, cycle 14). 
    C3 = 8 * CF;             % coupling strength in trun3 (maintenance, cycle 15). 
    C4 = 8 * CF;             % coupling strength in trun5-6 (larval eye disc). 
      
elseif plot_option == 5 || plot_option == 6 % eye 
    C1 = 0 * CF;             % coupling strength in trun1 (early cleavage cycles). 
    C2 = 0 * CF;             % coupling strength in trun2 (initiation, cycle 14). 
    C3 = 0 * CF;             % coupling strength in trun3 (maintenance, cycle 15). 
    C4 = 2.5 * CF;           % coupling strength in trun5-6 (larval eye disc). 

else error 'please enter plot_option between 1 and 7'
end


% The number and length of cell cycles are defined according to the cell type (somatic or germline)
% 1 iteration = 12 sec. 300 iterations = 1 hour. 
%--------------------------------------------------------------------------------------------------------

cell_type = 1;      % define the cell type. Choose 1 (somatic) or 2 (germline). Default = 1.
wipe = 0.5;         % proportion of nucleosomes that are converted to 'U' at replication. Default = 0.5

if cell_type == 1 % somatic
   
    teq = 50;           % number of iterations to pre- run for equilibration. Default = 50, this is cycle 1 for somatic.
    trun1 = 50;         % number of iterations per replication cycle in segment 1 (cycles 2-13). Default: 50 (= 10mins)
    trun2h = 300;       % segment 2h (use for early heatshock in cycle 14). Default: 300 (= 1 hour)
    trun2r = 500;       % segment 2r (the rest of cycle 14). Default: 500 (= 1h 40)
    trun2= trun2h+trun2r; % total number of iterations in segment 2. Cycle 14. Default: 800 (=2h 40)
    
    trun3h = 300;       % segment 3h (use for late heatshock in cycle 15.Default: 300 (= 1 hour)  
    trun3r = 500;        % segment 2r (the rest of cycle 15). Default: 500 (= 1h 40)
    trun3= trun3h+trun3r; % total number of iterations in segment 3. Cycle 15. Default: 800 (=2h 40)
    
 % currently the long simulations skip the rest of embryogenesis from 7.5h to 22h, starting straight in with 1st instar. To simulate the whole time course,
 % trun3r should be 4850 (total time up to trun3 = 22h)
    
    trun4h = 300;      % segment 4h  
    trun4r = 9*300;    % segment 4r. 4h and 4r are not used separately at the moment. Segment 4 is eye disc proliferation.
    trun4= trun4h+trun4r; % total number of iterations in segment 4. Default: 3000 (= 10 h). 

    trun5 = 8*300;   %8*300= 8 hours: ahead of furrow. Eya is activated in this window. At the end is the 2nd mitotic wave.
    trun6 =  20*300; %20*300= 20 hours: behind furrow. This brings you to the middle of third instar. 

    repeq = 1;        % number of replication cycles in equilibration phase. This is cycle 1 in somatic cells. Default = 1
    rep1 = 12;        % number of replication cycles in segment 1. These are the early cleavage cycles 2-12
    rep2 = 1;         % number of replication cycles in segment 2. Cycle 14
    rep3 = 1;         % number of replication cycles in segment 3. Cycle 15. Default = 1
    rep4 = 1;         % number of replication cycles in segment 4. This is the 1st instar larva, eye disc proliferation. Default = 7 
    rep5 = 1;         % mitosis after furrow (and mitotic wave)
    rep6 = 1;         % Postmitotic to end of 3rd instar there is no replication.

elseif cell_type == 2 % germline

    teq = 50;         % number of iterations for early cycles. This is the early cleavage cycles 1-9. Default = 50 (10 mins)
    trun1 = 200;       % number of iterations for cycle 10. Default: 200 (= 40mins)
    trun2h = 300;     % segment 2h. First part of cycle 11. Default: 300 (= 50 mins). Use for early heatshock, (2h10 to 3h10):
    trun2r = 50;     % segment 2r. Rest of cycle 11. Default: 50 (=10 mins)
    trun2= trun2h+trun2r; % total number of iterations in segment 2. Cycle 11. Default: 300 (= 1h)
    
    % total developmental time up to here: 9*50 + 500 = 950 iterations, 3h20mins. At
    % this stage the germ cells enter G2 arrest. There is one rep at 3h20
    % then no more!. 
    
    trun3h = 300;      % segment 3h. Start G1 arrest. default 300.
    trun3r = 1000+1150;      % segment 3r. Default 1000 
    trun3= trun3h+trun3r; % total number of iterations in segment 3. Cyc. Default: 1300 (= 4h20), if 7.5h total time course. 
    
    % currently the long simulations skip the rest of embryogenesis
    % from 7.5h to 22h.To simulate longer time course, extend trun3. 
    
    trun4h = 300;      % segment 4h  
    trun4r = 9*300;    % segment 4r. 4h and 4r are not used separately at the moment. Segment 4 is proliferation: can be used to represent mitosis during gonadogensis. Meiosis occurs later. 
    trun4= trun4h+trun4r; % total number of iterations in segment 4. Default: 3000 (= 10 h). 

    trun5 = 8*300;   %8*300= 8 hours: 
    trun6 =  20*300; %20*300= 20 hours: This brings you to the middle of third instar. 

    repeq = 9;       %early cleavage.
    rep1 = 1;        % number of replication cycles in segment 1. Cyle 10
    rep2 = 1;        % number of replication cycles in segment 2. Cycle 11.
    rep3 = 1;        % number of replication cycles in segment 3. Cycle 12. G2 arrest.
    rep4 = 1;        % number of replication cycles in segment 4. This is the 1st instar larva, mitotic proliferation. Default = 7
    rep5 = 1;        % mitosis after furrow (and mitotic wave)
    rep6 = 1;        % Postmitotic to end of 3rd instar, there is no replication.
    
else error 'please enter 1 (somatic) or 2 (germline) for cell_type'
        
end

if plot_option <= 4 || plot_option == 7 % short runs
    rep = repeq+rep1+rep2+rep3; % this goes only up to the end of embryogenesis. 
    t = teq*repeq+ trun1*rep1 + trun2*rep2 +trun3*rep3;

elseif plot_option == 5 || plot_option == 6 % long runs
    rep = repeq+rep1+rep2+rep3+rep4+rep5+rep6; % this runs the larva up to end of third instar
    t = teq*repeq + trun1*rep1 + trun2*rep2 + trun3*rep3 + trun4*rep4 + trun5*rep5 + trun6 *rep6;
    
% these are the inputs for the Hill function in section 5 (up) and 6 (down).
    
    p1up = 0.6;        
    p1down = 0.25; 
    C_up = 0.4;
    h_up = 3; 
    C_down = 0.85;
    h_down = 2;        
 

 else error 'please enter plot_option 1 to 7'
 end
    
% Set up the inputs for the promoter (p1 (TF binding) and p2 (TF dissociation)
%----------------------------------------------------------------------------
% p2 does not change, but it could be set to different values, currently the same for all plot options
%------------------------------------------------------------------------------------------------------------
p2Val =  [0.1, 0.1];  % p2 values, segment 1. 
p2hVal = p2Val;    % p2 values for heatshock in segment 2. 
p2aVal = p2Val;    % p2 values, segment 2. 
p2h2Val = p2Val;    % p2 values for heatshock in segment 3. 
p2bVal = p2Val;      % p2 values, segment 3. 
p2h3Val = p2Val;    % p2 values for heatshock in segment 4. 
p2cVal = p2Val;      % p2 values, segment 4. 
p2dVal = p2Val;      % p2 values, segment 5. 
p2eVal = p2Val;      % p2 values, segment 5. 

% p1 is different for different plot options.  
%-------------------------------------------------------------------------------------

 if plot_option == 1 || plot_option == 2 % simulates promoter activity in anterior (A) and posterior (P) over time. 

    p1Val =  [0.001, 0.001];    % p1 values, segment 1. same for all. Default [0.001, 0.001]
    p1hVal = [0.001, 0.6];      % p1 values for cycle 14 [A, P]. Default [0.001, 0.6]
    p1aVal = [0.001, 0.6];      % p1 values for  cycle 14 [A, P]. Default [0.001, 0.6]
    p1h2Val = [0.6, 0.6];       % p1 values for cycle 15 [A, P]. Default [0.6, 0.6]
    p1bVal = [0.6, 0.6];        % p1 values for cycle 15 [A, P]. Default [0.6, 0.6]
    p1h3Val = [0.6, 0.6];       % p1 values for heatshock in segment4. Not used in this time course but needs input, value not important  
    p1cVal = [0.6, 0.6];        % p1 values, segment 4. Not used in this time course but needs input, value not important 
    p1dVal = [0.6, 0.6];        % p1 values, segment 5. Not used in this time course but needs input, value not important 
    p1eVal = [0.6, 0.6];        % p1 values, segment 6. Not used in this time course but needs input, value not important     
    
 elseif plot_option == 3 || plot_option == 4 || plot_option == 7 % simulates heat shock in both compartments over time. Option for no, early or late heatshock.

    p1Val =  [0.001, 0.001];    % p1 values, segment 1. same for all. Default[0.001, 0.001] 
    p1hVal = [0.06, 0.06];      % p1 values for cycle 14. No heatshock, basal promoter activity. Default[0.06, 0.06] 
    p1aVal = [0.06, 0.06];      % p1 values for cycle 14    " 
    p1h2Val = [0.06, 0.06];     % p1 values for cycle 15    "
    p1bVal = [0.06, 0.06];      % p1 values for cycle 15    "
    p1h3Val = [0.06, 0.06];     % p1 values for heatshock in segment4.Not used in this time course but needs input, value not important  
    p1cVal = [0.06, 0.06];      % p1 values, segment 4. Not used in this time course but needs input, value not important 
    p1dVal = [0.06, 0.06];      % p1 values, segment 5. Not used in this time course but needs input, value not important 
    p1eVal = [0.06, 0.06];      % p1 values, segment 6. Not used in this time course but needs input, value not important 
    
        if heatshock == 0
            p1hVal = [0.06, 0.06];  % no heatshock
        elseif heatshock == 1
            p1hVal = [0.6, 0.6];    % early heatshock
        elseif heatshock == 2
            p1h2Val = [0.6, 0.6];   % late heatshock
        else
            error 'please enter heatshock = 0(none), 1(early) or 2(late)'
        end
    
elseif plot_option == 5 || plot_option == 6 % eye disc, eya PRE
    p1Val =  [0.001, 0.001];    % p1 values, segment 1. 
    p1hVal = [0.01, 0.01];      % p1 values for cycle 14. specific for 1,2 and 3,4 and 5,6 (and 7)
    p1aVal = [0.01, 0.01];      % p1 values, cycle 14. specific for 1,2 and 3,4 and 5,6 (and 7)
    p1h2Val = [0.01, 0.01];       % p1 values for cycle 15 [A, P]
    p1bVal = [0.01, 0.01];        % p1 values, cycle 15 [A,P
    p1h3Val = [0.01, 0.01];      % p1 values for heatshock in segment4. 
    p1cVal = [0.01, 0.01];     % p1 values, segment 4. 
    p1dVal =  [p1up,p1up];%[0.32, 0.32];     % p1 values, segment 5. Use eya values if running plot option 7 or 8.
    p1eVal =  [p1down,p1down];%[0.3, 0.3];     % p1 values, segment 6. Use eya values if running plot option 7 or 8.
 
else
    error 'please enter plot_option 1-7'

end

V = numel (p1Val);
  
% PRE/TRE inputs (p3, p4, p5, p6, m, tf) 
%-------------------------------------------------------------------------------------
if plot_option == 1 || plot_option == 3  % bxd PRE
    p5 = 0.04; % 0.04 
elseif plot_option == 5 || plot_option == 7 %eya PRE or bxd PRE.
    p5 = p5in;
end
         
if plot_option == 1 || plot_option == 3 || plot_option == 5 || plot_option == 7 % simple time course
    p3ser = [0.25,0.25]; % default [0.25,0.25], A and P.
    p4ser = p3ser;     
        
    if PcG == 1
        p3_3ser = p3ser; % values during the maintenance phase: change p3_3ser to [0.001, 0.001]for PcG mutant.
    elseif PcG == 0
        p3_3ser = [0.001,0.001]; % values during the maintenance phase: change p3_3ser to [0.001, 0.001]for PcG mutant.
    end
        
    if TrxG == 1
        p4_3ser = p4ser; % values during the maintenance phase: change p4_3ser to [0.001, 0.001]for TrxG mutant.
    elseif TrxG == 0
        p4_3ser = [0.001,0.001]; % values during the maintenance phase: change p3_3ser to [0.001, 0.001]for PcG mutant.
    end
        
    p5ser = [p5,p5];  
    p6 = 0.5;           % probability of going to A or M from U in noisy conversion. Default = 0.5                                  
    mser = [m, m];  
    tfser = [tf,tf]; 

 elseif plot_option == 2 || plot_option == 4 % parameter space memory of silencing (2) or activation (4) 
     p3ser = [0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,0.9]; % default [0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,0.9];
     p4ser = p3ser;     
     p3_3ser = p3ser;
     p4_3ser = p4ser;
     p5ser = [0.00125,0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.32]; % default [0.00125,0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.32];
     p6 = 0.5;                                     
     mser = [m,m,m,m,m,m,m,m,m];  
     tfser = [tf, tf, tf, tf, tf, tf, tf, tf,tf]; 

 elseif plot_option == 6 % eye disc parameter space series.    
     p3ser = [0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.4,0.6]; % default [0.005,0.01,0.05,0.1,0.15,0.25,0.3,0.4,0.8];
     p4ser = p3ser;    
     p3_3ser = p3ser;
     p4_3ser = p4ser;
     p5ser = [0.04,0.06,0.08,0.1,0.12,0.17,0.2,0.4,0.8];   % default [0.04,0.06,0.08,0.1,0.12,0.16,0.2,0.4,0.8];
     p6 = 0.5;                                     
     mser = [m,m,m,m,m,m,m,m,m];  % no. of nucleosomes in PRE
     tfser = [tf, tf, tf, tf, tf, tf, tf, tf,tf]; % no. of sites in enhancer

 else error 'please enter plot_option 1 to 7'
  
 end
  
parameter1 =  p3ser;    % assign parameter 1 to the series that is varied first.
parameter2 =  p5ser;    % assign parameter 2 to the series that is varied second.

% Define coupling inputs according to plot option.
%-------------------------------------------------------------------------------------
if plot_option == 1 || plot_option == 3 || plot_option == 7 
        Cp_eser1 = [C1,C1];       % PRE-enhancer, segment 1. 
        Ce_pser1 = Cp_eser1;    % enhancer-PRE, segment 1.
        Cp_eser2 = [C1,C2];      % PRE-enhancer, segment 2.
        Ce_pser2 = Cp_eser2;    % enhancer-PRE, segment 2.         
        Cp_eser3 = [C1,C3];      % PRE-enhancer, segment 3. 
        Ce_pser3 = Cp_eser3;    % enhancer-PRE, segment 3.

 elseif plot_option == 2 || plot_option == 4 
        Cp_eser1 =  [C1,C1,C1,C1,C1,C1,C1,C1,C1];   % PRE-enhancer, segment 1. 
        Ce_pser1 = Cp_eser1;               % enhancer-PRE, segment 1.
        Cp_eser2 = [0,C2,C2,C2,C2,C2,C2,C2,C2];    % PRE-enhancer, segment 2. The 0 in first position is used to generate controls
        Ce_pser2 = Cp_eser2;               % enhancer-PRE, segment 2.         
        Cp_eser3 = [0,C3,C3,C3,C3,C3,C3,C3,C3];    % PRE-enhancer, segment 3. The 0 in first position is used to generate controls 
        Ce_pser3 = Cp_eser3;               % enhancer-PRE, segment 3.
    
 elseif plot_option == 5 
        Cp_eser1 =  [C1,C1];   % PRE-enhancer, segment 1. 
        Ce_pser1 = Cp_eser1;               % enhancer-PRE, segment 1.
        Cp_eser2 = [C1,C1];    % PRE-enhancer, segment 2.
        Ce_pser2 = Cp_eser2;               % enhancer-PRE, segment 2.         
        Cp_eser3 = [C1,C4];    % PRE-enhancer, segment 3. 
        Ce_pser3 = Cp_eser3;               % enhancer-PRE, segment 3.

 elseif plot_option == 6
        Cp_eser1 =  [C1,C1,C1,C1,C1,C1,C1,C1,C1];   % PRE-enhancer, segment 1. 
        Ce_pser1 = Cp_eser1;               % enhancer-PRE, segment 1.
        Cp_eser2 = [C1,C1,C1,C1,C1,C1,C1,C1,C1];    % PRE-enhancer, segment 2.
        Ce_pser2 = Cp_eser2;               % enhancer-PRE, segment 2.         
        Cp_eser3 = [C4,C4,C4,C4,C4,C4,C4,C4,C4];    % PRE-enhancer, segment 3. 
        Ce_pser3 = Cp_eser3;               % enhancer-PRE, segment 3.

 else
        error 'please enter plot_option 1-7'
 end
    
% Generate the input vectors from the parameter series "ser" using function "memory_inputs"
%------------------------------------------------------------------------------------
 if plot_option == 2 || plot_option == 4 || plot_option == 6 %parameter space versions
        [p3input,p4input,p3_3input,p4_3input,minput,p5input,Cp_einput1,Ce_pinput1,Cp_einput2,Ce_pinput2,Cp_einput3,Ce_pinput3,tfinput,parameter2input] = memory_inputs(parameter1,parameter2,p3ser,p4ser,p3_3ser,p4_3ser,mser,p5ser,Cp_eser1,Ce_pser1,Cp_eser2,Ce_pser2,Cp_eser3,Ce_pser3,tfser);
        P3Val = p3input; 
        P4Val = p4input;
        P3_3Val = p3_3input; 
        P4_3Val = p4_3input;
        mVal = minput;
        p5Val= p5input;

        Cp_eVal1 = Cp_einput1;
        Ce_pVal1 = Ce_pinput1;
        Cp_eVal2 = Cp_einput2;
        Ce_pVal2 = Ce_pinput2;
        Cp_eVal3 = Cp_einput3;
        Ce_pVal3 = Ce_pinput3;
        tfVal = tfinput;

 else % all other plot options run a simple time course and take the inputs directly from the entered series.
        P3Val = p3ser; 
        P4Val = p4ser;
        P3_3Val = p3_3ser; 
        P4_3Val = p4_3ser;
        mVal = mser;
        p5Val= p5ser;

        Cp_eVal1 = Cp_eser1;
        Ce_pVal1 = Ce_pser1;
        Cp_eVal2 = Cp_eser2;
        Ce_pVal2 = Ce_pser2;
        Cp_eVal3 = Cp_eser3;
        Ce_pVal3 = Ce_pser3;
        tfVal = tfser;    
 end

C = numel(P3Val);

% Calculate window size for evaluation of means and lifetimes
%------------------------------------------------------------------------------------------------------
     
windowsize = 50;%  the window calculates one mean per 10 mins.
step = windowsize;%step to move window 
wpositions = (1:step:t);
W = numel (wpositions);
    
numstepeq = teq*repeq/windowsize;
numstep1 = trun1*rep1/windowsize;
numstep2 = trun2*rep2/windowsize;
numstep2h = trun2h*rep2/windowsize;
numstep3 = trun3*rep3/windowsize;
numstep4 = trun4*rep4/windowsize;
numstep5 = trun5*rep5/windowsize;
numstep6 = trun6*rep6/windowsize;

         
if plot_option <= 4 || plot_option == 7 %the short time courses
        numstep = numstepeq+numstep1+numstep2+numstep3;

elseif plot_option == 5 || plot_option == 6 % the long time courses
        numstep = numstepeq+numstep1+numstep2+numstep3+numstep4+numstep5+numstep6;    
   
end
    
numstephs = numstepeq+numstep1+numstep2h; %number of windows to end of early heatshock.
numstepAfterhs = numstep - numstephs; %number of windows after heatshock to end of time course
     
    
% set up output matrices for the repetitions of whole simulatin (k loop).
%---------------------------------------------------------------------------------------------------------- 
        
sum_timeA = zeros (t, V*C);
sum_timeM = zeros (t, V*C);
sum_timeB = zeros (t, V*C);
sum_timeF = zeros (t, V*C);
    
ALL_stateP = zeros (W, 2*n); % States for each window.
ALL_stateE = zeros (W, 2*n);
    
ALL_stateP_no = zeros (W, n);
ALL_stateP_with = zeros (W, n);
ALL_stateE_no = zeros (W, n);
ALL_stateE_with = zeros (W, n);
    
ALL_life = zeros (6, 2*n); % Lifetime of each state, in cell cycles. Rows 1-6 = A,U,M,B,N,F. 
       
BIG_Cp = zeros (t*n,4*V*C); % this has ALL the coupling parameters, including teq
BIG_Pco = zeros (t*n,4*V*C); % this has ALL the parameters p1co to p4 co
        
BIG_timeA = zeros (t*n, V*C); %This has ALL the data for the time courses, INCLUDING the teq.
BIG_timeM = zeros (t*n, V*C);
BIG_timeB = zeros (t*n, V*C);
BIG_timeF = zeros (t*n, V*C);
    
BIG_meanA = zeros (W*n, V*C); %This has ALL the data for the means including teq.
BIG_meanM = zeros (W*n, V*C);
BIG_meanB = zeros (W*n, V*C);
BIG_meanF = zeros (W*n, V*C);

for k = 1:n

    % Set up output matrices for the coupling series (f loop). 
    %------------------------------------------------------------------------------------------------------ 
    
    timeA_ALL = zeros (t, V*C);
    timeM_ALL = zeros (t, V*C);
    timeB_ALL = zeros (t, V*C);
    timeF_ALL = zeros (t, V*C);

    meanA_ALL = zeros (W, V*C);
    meanM_ALL = zeros (W, V*C);
    meanB_ALL = zeros (W, V*C);
    meanF_ALL = zeros (W, V*C);
    
    stateP_ALL = zeros (W, V*C);
    stateE_ALL = zeros (W, V*C);
    
    cyc_lifetimes_ALL = zeros (6, V*C);   
        
    Cp_ALL = zeros (t,4*V*C);
    Pco_ALL = zeros (t,4*V*C);

for f = 1:C % runs the j block C times. For time courses, C = 2 (with and without coupling). For parameter space plots, C = 81 (9*9 values of parameters 1 and 2).
       
    p3 = P3Val(f); % constant in time courses, varied in parameter space search (parameter 1)
    p4 = P4Val(f); % "
    
    p3_3 = P3_3Val(f);
    p4_3 = P4_3Val(f);
    
    
    m  = mVal(f);  % constant in all
    p5 = p5Val(f); %constant in time courses, varied in parameter space search (parameter 2)
    
    Cp_e1 = Cp_eVal1(f); %constant in parameter space search, varied in time courses. (0 or value).
    Ce_p1 = Ce_pVal1(f); % "
    
    Cp_e2 = Cp_eVal2(f); % "
    Ce_p2 = Ce_pVal2(f); % "
    
    Cp_e3 = Cp_eVal3(f); % "
    Ce_p3 = Ce_pVal3(f); % "
    
    tf  = tfVal(f); % constant in all
       
        
% set up output matrices for the p1 series (j loop).
%----------------------------------------------------------------------------------------------------------

    timeA = zeros (t,V);  
    timeM = zeros (t,V);  
    timeB = zeros (t,V);  
    timeF = zeros (t,V);
    
    cyc_meanA = zeros (W, V);
    cyc_meanM = zeros (W, V);
    cyc_meanB = zeros (W, V);
    cyc_meanF = zeros (W, V);
    
    cyc_stateP = zeros (W, V);
    cyc_stateE = zeros (W, V);

    cyc_lifetimes = zeros (6, V);
    
    Cp = zeros (t,4*V);
    Pco = zeros (t,4*V);
    
    
for j = 1:V % runs the time course V times, for different input p1 values (for plot option 1 and 2, this is anterior, posterior)
    
       p1 = p1Val(j); % input for "timecourse_function" during first replication cycles 
       p1h = p1hVal(j); % input for "timecourse_function" during first heatshock 
       p1a = p1aVal(j);% input for timecourse_function during second replication cycles 
       p1h2 = p1h2Val(j); % input for "timecourse_function" during second heatshock 
       p1b = p1bVal(j);% input for timecourse_function during third replication cycles
       p1h3 = p1h3Val(j); % input for "timecourse_function" during third heatshock 
       p1c = p1cVal(j);% input for timecourse_function during fourth replication cycles
       p1d = p1dVal(j);% input for timecourse_function during fifth replication cycles
       p1e = p1eVal(j);% input for timecourse_function during sixth replication cycles

       p2 = p2Val(j); % input for "timecourse_function" during first replication cycles 
       p2h = p2hVal(j); % input for "timecourse_function" during first heatshock
       p2a = p2aVal(j);% input for timecourse_function during second replication cycles 
       p2h2 = p2h2Val(j); % input for "timecourse_function" during second heatshock 
       p2b = p2bVal(j);% input for timecourse_function during third replication cycles
       p2h3 = p2h3Val(j); % input for "timecourse_function" during third heatshock
       p2c = p2cVal(j);% input for timecourse_function during fourth replication cycles
       p2d = p2dVal(j);% input for timecourse_function during fifth replication cycles
       p2e = p2eVal(j);% input for timecourse_function during sixth replication cycles
 
% run the simulation to equilibrate the PRE and the enhancer with "timecourse_function"
%---------------------------------------------------------------------------------------       




    if PRE_init == ('M')%defines the first row of nucs according to the input for the initial state.
        nucstart = 1; 
    elseif PRE_init == ('A')
        nucstart = -1;
    elseif PRE_init == ('U')
        nucstart = 0;
    else error 'PRE initial state must be M, A or U'
    end

    if Enhancer_init == ('B')%defines the first row of sites according to the input for the initial state.
        sitestart=1; 
    elseif Enhancer_init == ('F')
        sitestart = 0;
    else error 'Enhancer initial state must be B or F'
    end
    
    
couple_ep_p3start= 1; %set to 1 for the first run, then is calculated from the enhancer state.
couple_ep_p4start= 1; %set to 1 for the first run, then is calculated from the enhancer state.

Cp_e = Cp_e1;
Ce_p = Ce_p1;

    for ee = 1:repeq
        
        ti=teq;

        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3,p4,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start); 
        [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);%these are given to the next timecourse function

        hi = (ee*(ti-1))+ee;
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTeq(lo:hi,:)= OUT;% puts OUT into OUTeq1 at given position.  
    
    end
%OUTeq = OUT;
% set inputs for the next run (this is the real start of the simulation
% after equilibration. Takes the nucs and sites from the equilibration
% run.)

% run the simulation with replication for the PRE and the enhancer with "timecourse_function"
%---------------------------------------------------------------------------------------------    

%nucstart= nucs(teq,:); %defines the first row of nucs as the last row of the last one. 
%sitestart= sites(teq,:);%defines the first row of sites

couple_ep_p3start= OUTeq (teq*repeq,13); %last value calculated from the enhancer state.
couple_ep_p4start= OUTeq (teq*repeq,14); %last value calculated from the enhancer state.

Cp_e = Cp_e1;
Ce_p = Ce_p1;

OUTrep1 = zeros (trun1*rep1,18); % set up output matrix "OUTrep1"

    for e = 1:rep1

        ti=trun1;
   
        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3,p4,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
        [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);%these are given to the next timecourse function

    
        hi = (e*(ti-1))+e;
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTrep1(lo:hi,:)= OUT;% puts OUT into OUTrep1 at given position.   
    
    end
    
    couple_ep_p3start= OUTrep1 (ti*rep1,13); %last value calculated from the enhancer state.
    couple_ep_p4start= OUTrep1 (ti*rep1,14); %last value calculated from the enhancer state. 
    
    Cp_e = Cp_e2;
    Ce_p = Ce_p2;
    
    
    
    OUTrep2 = zeros (trun2*rep2,18); % set up output matrix "OUTrep2"
    
    for d = 1:rep2

        ti=trun2h;
    
        p1 = p1h;%change the p1 value;
        p2 = p2h;%change the p1 value;
    
   
        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3,p4,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
        hi = ti + (trun2*(d-1));
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTrep2(lo:hi,:)= OUT;% puts OUT into OUTreph at given position.  
    
        nucstart= nucs (ti,:);  
        sitestart = sites (ti,:);
        couple_ep_p3start= OUTrep2 (ti,13); %last value calculated from the enhancer state.
        couple_ep_p4start= OUTrep2 (ti,14); %last value calculated from the enhancer state.
    
        ti=trun2r;
    
        p1 = p1a;%change the p1 value;
        p2 = p2a;%change the p1 value;
    
        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3,p4,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
        [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);
 
        hi = trun2 + (trun2 * (d-1));
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTrep2(lo:hi,:)= OUT;% puts OUT into OUTrep2 at given position.
     
    end
    
    couple_ep_p3start= OUTrep2 (trun2*rep2,13); %last value calculated from the enhancer state.
    couple_ep_p4start= OUTrep2 (trun2*rep2,14); %last value calculated from the enhancer state. 
    
    Cp_e = Cp_e2;
    Ce_p = Ce_p2;
   
    
    OUTrep3 = zeros (trun3*rep3,18); % set up output matrix "OUTrep3"
    
    for c = 1:rep3 
   
        ti=trun3h;
    
        p1 = p1h2;%change the p1 value;
        p2 = p2h2;%change the p1 value;       
        
        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
        hi = ti + (trun3* (c-1));
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTrep3(lo:hi,:)= OUT;% puts OUT into OUTrep3 at given position.  
    
        nucstart= nucs (ti,:);  
        sitestart = sites (ti,:);
        couple_ep_p3start= OUTrep3 (ti,13); %last value calculated from the enhancer state.
        couple_ep_p4start= OUTrep3 (ti,14); %last value calculated from the enhancer state.
    
        ti=trun3r;
    
        p1 = p1b;%change the p1 value;
        p2 = p2b;%change the p1 value;
                    
        
        [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
        [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);
   
        hi = trun3 + (trun3 * (c-1));
        lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
        OUTrep3(lo:hi,:)= OUT;% puts OUT into OUTrep1 at given position.
    
    end
    
    
    if plot_option <= 4 || plot_option == 7 % the short runs.
        
        OUTALL=[OUTeq;OUTrep1;OUTrep2;OUTrep3];
    
    elseif plot_option == 5 || plot_option == 6 % the long runs
            
        couple_ep_p3start= OUTrep3 (trun3*rep3,13); %last value calculated from the enhancer state.
        couple_ep_p4start= OUTrep3 (trun3*rep3,14); %last value calculated from the enhancer state. 
    
        Cp_e = Cp_e2;
        Ce_p = Ce_p2;
         
        
        OUTrep4 = zeros (trun4*rep4,18); % set up output matrix "OUTrep3"
            
    for b = 1:rep4

            ti=trun4h;
    
            p1 = p1h3;%change the p1 value;
            p2 = p2h3;%change the p1 value;
    
   
            [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
            hi = ti + (trun4* (b-1));
            lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
   
            OUTrep4(lo:hi,:)= OUT;% puts OUT into OUTreph at given position.  
    
            nucstart= nucs (ti,:);  
            sitestart = sites (ti,:);
            
            couple_ep_p3start= OUTrep4 (ti,13); %last value calculated from the enhancer state.
            couple_ep_p4start= OUTrep4 (ti,14); %last value calculated from the enhancer state.
    
            ti=trun4r;
    
            p1 = p1c;%change the p1 value;
            p2 = p2c;%change the p1 value;
    
    
            [nucs, sites, OUT] = timecourse_function(ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
            [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);
 
            hi = trun4 + (trun4 * (b-1));
            lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
  
            OUTrep4(lo:hi,:)= OUT;% puts OUT into OUTrep2 at given position
    end
            
        nucstart = nucs (ti-1,:); %ignores the last rep values of run4
        sitestart = sites (ti-1,:);
        
        couple_ep_p3start= OUTrep4 (trun4*rep4,13); %last value calculated from the enhancer state.
        couple_ep_p4start= OUTrep4 (trun4*rep4,14); %last value calculated from the enhancer state. 
    
        Cp_e = Cp_e3;
        Ce_p = Ce_p3;
         
        
        OUTrep5 = zeros (trun5*rep5,18); % set up output matrix "OUTrep5"
            
    for a = 1:rep5

            ti=trun5;
    
            p1 = p1d;%change the p1 value;
            p2 = p2d;%change the p1 value;
            
    
            [nucs, sites, OUT,UPHILLOUT] = timecourse_functionUPHILL(C_up,h_up,ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
            [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);
 
            hi = trun5 + (trun5 * (a-1));
            lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
  
            OUTrep5(lo:hi,:)= OUT;% puts OUT into OUTrep5 at given position
            %end of new part.
            
   end
            
        couple_ep_p3start= OUTrep5 (trun5*rep5,13); %last value calculated from the enhancer state.
        couple_ep_p4start= OUTrep5 (trun5*rep5,14); %last value calculated from the enhancer state. 
    
        Cp_e = Cp_e3;% can also change this to a version 4.
        Ce_p = Ce_p3;
         
        
        OUTrep6 = zeros (trun6*rep6,18); % set up output matrix "OUTrep5"
            
   for aa = 1:rep6

            ti=trun6;
    
            p1 = p1e;%change the p1 value;
            p2 = p2e;%change the p2 value;
    
            [nucs, sites, OUT, DOWNHILLOUT] = timecourse_functionDOWNHILL(C_down,h_down,ti,m,tf,model,p3_3,p4_3,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start);
    
            [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe);
 
            hi = trun6 + (trun6 * (aa-1));
            lo = (hi-(ti-1));% calculates lower and upper values for intervals of output matrix.
  
            OUTrep6(lo:hi,:)= OUT;% puts OUT into OUTrep5 at given position
            %end of new part.
            
    end
            OUTALL=[OUTeq;OUTrep1;OUTrep2;OUTrep3;OUTrep4;OUTrep5;OUTrep6];
        
        
    else error 'please enter plot_option 1-7'
        
    end
                   
    %put result of each run in timeX.
    %------------------------------------------------------------------------------------------

    timeA (:,j) = OUTALL(1:t,5);  %timecourse A: This is the raw data for each iteration. (number of 'A' /m)
    timeM (:,j) = OUTALL(1:t,4);  %timecourse M
    timeB (:,j) = OUTALL(1:t,9);  %timecourse B
    timeF (:,j) = OUTALL(1:t,10); %timecourse F
    
    %calculate means and states for each run with "statefunction". (similar to plot 1 and 2). 
    %------------------------------------------------------------------------------------------

    Arun = OUTALL (1:t,5);
    Mrun = OUTALL (1:t,4);
    Brun = OUTALL (1:t,9);
    Frun = OUTALL (1:t,10);
 
    [WA, WM, WB, WF, stateP, stateE] = statefunction (windowsize,t,step,W,Arun,Mrun,Brun,Frun);

    cyc_meanA (:,j)= WA(:,1); % stores the mean proportion of "A" states across a time window (given by 'windowsize')
    cyc_meanM (:,j)= WM(:,1); % same for M
    cyc_meanB (:,j)= WB(:,1); % same for B
    cyc_meanF (:,j)= WF(:,1); % same for F

    cyc_stateP (:,j)= stateP(:,1); % stores the mean state of the system each time window for the PRE.
    cyc_stateE (:,j)= stateE(:,1); % stores the mean state of the system each time window for the PRE.

%Calculate the counts and transitions from stateP and stateE using "lifetimefunction". 
%This will give lifetimes in terms of windows: for how many time windows does a given state survive. 
%-----------------------------------------------------------------------------------------------

    state = stateP; [counts] = lifetimefunction(state, W);
 
    lifeA = (counts (W,1))/(counts (W,2));
    lifeM = (counts (W,3))/(counts (W,4));
    lifeU = (counts (W,5))/(counts (W,6));
 
    state = stateE; [counts] = lifetimefunction(state, W);
 
    lifeB = (counts (W,1))/(counts (W,2));
    lifeF = (counts (W,3))/(counts (W,4));
    lifeN = (counts (W,5))/(counts (W,6));
 

% put the outputs into the matrix "cyc_lifetimes" for each iteration of j.
%----------------------------------------------------------------------------
 
    cyc_lifetimes (1,j)= lifeA;
    cyc_lifetimes (2,j)= lifeU;
    cyc_lifetimes (3,j)= lifeM;
    cyc_lifetimes (4,j)= lifeB;
    cyc_lifetimes (5,j)= lifeN;
    cyc_lifetimes (6,j)= lifeF;

    cyc_lifetimes(isnan(cyc_lifetimes)) = 0;
    
   
    %Store the parameters in blocks of 4, one for each iteration of j.(couple
    %p1-p4, p1-p4co.)
    %------------------------------------------------------------------------------------------
    
    Cp (:,j*4-3:j*4)= OUTALL (:,11:14);%couple_p1; couple_p2; couple_ep_p3;couple_ep_p4;
    Pco (:,j*4-3:j*4)= OUTALL (:,15:18);%p1co, p2co, p3co, p4co
    

end % end j loop (p1 values- 1:V).
   
    block = (V*f)-V+1;
    
    timeA_ALL (:,block:V*f) = timeA (:,:); 
    timeM_ALL (:,block:V*f) = timeM (:,:);
    timeB_ALL (:,block:V*f) = timeB (:,:);
    timeF_ALL (:,block:V*f) = timeF (:,:);

    meanA_ALL (:,block:V*f) = cyc_meanA (:,:);
    meanM_ALL (:,block:V*f) = cyc_meanM (:,:);
    meanB_ALL (:,block:V*f) = cyc_meanB (:,:);
    meanF_ALL (:,block:V*f) = cyc_meanF (:,:);

    stateP_ALL (:,block:V*f) = cyc_stateP (:,:);
    stateE_ALL (:,block:V*f) = cyc_stateE (:,:);

    cyc_lifetimes_ALL (:,block:V*f)= cyc_lifetimes (:,:);
    
    block1 = (4*V*f)-(4*V-1); 
    
    Cp_ALL (:,block1:4*V*f)= Cp (:,:); %puts the parameters in sets of four for each p1 value, consecutive blocks of p1 values.
    Pco_ALL (:,block1:4*V*f) = Pco (:,:);   
   
end %end f loop (coupling values and p3 p4)
    
    % sum up all repetitions for timecourse. The mean is calculated in the
    % function 'time_state'
    
    sum_timeA = sum_timeA + timeA_ALL;
    sum_timeM = sum_timeM + timeM_ALL;
    sum_timeB = sum_timeB + timeB_ALL;
    sum_timeF = sum_timeF + timeF_ALL;
    
    % extract state information (used in time_state). 
    
    ALL_stateP_no(:,k) = stateP_ALL(:,1); ALL_stateP_with(:,k) = stateP_ALL(:,3); %makes 2 separate files, with and without coupling.ONLY FOR THE FIRST P1 VALUE!!
    ALL_stateE_no(:,k) = stateE_ALL(:,1); ALL_stateE_with(:,k) = stateE_ALL(:,3); % makes 2 separate files, with and without coupling. ONLY FOR THE FIRST P1 VALUE!!
    
    % extract lifetime information 
  
    ALL_life (:,2*k-1)  = cyc_lifetimes_ALL(:,1);ALL_life (:,2*k) = cyc_lifetimes_ALL(:,3);%this is not used in any plot at the moment.
       
    % add all timecourse outputs vertically 
  
    hit = k*t;
    lot = hit +1 - t;
        
    BIG_Cp (lot:hit,:) = Cp_ALL (:,:); %add all coupling parameters to matrix vertically.
    BIG_Pco (lot:hit,:) = Pco_ALL (:,:); %add all p1co-p4co to matrix vertically.
 
    BIG_timeA (lot:hit,:) = timeA_ALL (:,:); %add all timecourses to matrix vertically.
    BIG_timeM (lot:hit,:) = timeM_ALL (:,:);
    BIG_timeB (lot:hit,:) = timeB_ALL (:,:);
    BIG_timeF (lot:hit,:) = timeF_ALL (:,:);
    
    him = k*W;
    lom = him+1-W;
    
    BIG_meanA (lom:him,:) = meanA_ALL (:,:); %add them all up vertically
    BIG_meanM (lom:him,:) = meanM_ALL (:,:); %add them all up 
    BIG_meanB (lom:him,:) = meanB_ALL (:,:); %add them all up 
    BIG_meanF (lom:him,:) = meanF_ALL (:,:); %add them all up 
    
    BIG_stateP (lom:him,:) = stateP_ALL (:,:); %add them all up vertically
    BIG_stateE (lom:him,:) = stateE_ALL (:,:); %add them all up 


end %k loop. (k=1:n) repetitions for mean


% Call plot functions. 
% -----------------------------------------------------------------------

    
if plot_option == 1||plot_option == 3 % memory of silencing or activation time course. 
        
   [Bend, Aend, AMend,meanBend, stdBend, meanAend,stdAend,meanAMend,stdAMend,...
     memory_score_nw, mean_memory_nw, std_memory_nw]...
    = plot_option_1_3(note,plot_option,BIG_meanA,BIG_meanM,BIG_meanB,n,numstepeq,numstep1,numstep2,numstep3,numstep);
        
elseif plot_option == 2 % memory of silencing parameter space.
        
   [One_diffBm2,meanLRBm,LRmeanBm0,singlevalBm]...
       = plot_option_2(parameter1,parameter2,p3ser,p3_3ser,p4_3ser,p4ser,p5ser,mser,tfser,...
       Cp_eser1,Cp_eser2,Cp_eser3,Ce_pser2,Ce_pser3,n,...
       V,C,BIG_meanB,BIG_meanA,numstep,numstep3,note);
  
           
elseif plot_option == 4 %memory of activation parameter space.
        
   [One_diffBm2,meanLRBm,LRmeanBm0,singlevalBm]...
       = plot_option_4(parameter1,parameter2,p3ser,p3_3ser,p4_3ser,p4ser,p5ser,mser,tfser,...
       Cp_eser1,Cp_eser2,Cp_eser3,Ce_pser2,Ce_pser3,n,...
       V,C,BIG_meanB,BIG_meanA,numstep,numstepAfterhs,note);
    
elseif plot_option == 5 %eye disc plots.
        
   [Timeno, Timewith, TimewithS,Discno, Discwith, AMDiscno, AMDiscwith, STDno, STDwith]...
       = plot_option_5(BIG_meanA,BIG_meanM,BIG_meanB,n,...
       numstepeq,numstep1,numstep2,numstep3,numstep4,numstep5,numstep6,numstep, note, PRE_ID);

elseif plot_option == 6  %eye disc parameter space.
            
   [Index, data_with,MeanBdisc,MeanAMdisc, P1APmean, FitP1A_All, FitP1P_All,Varlist]...
       = plot_option_6(parameter1,parameter2,p3ser,p3_3ser,p4_3ser,p4ser,p5ser,mser,tfser,...
       Cp_eser1,Cp_eser2,Cp_eser3,Ce_pser2,Ce_pser3,P3Val,p5Val,n,...
       BIG_meanA,BIG_meanM,BIG_meanB,numstepeq,numstep1,numstep2,numstep3,numstep4,numstep5,numstep6,numstep,note);

elseif plot_option == 7 % Time course, single runs and means; state analysis for last 10 mins.
        
   [figure1] = plot_option_7(t,k,n,sum_timeA,sum_timeB,sum_timeM,sum_timeF,...
       ALL_stateP_no,ALL_stateP_with,ALL_stateE_no,ALL_stateE_with,note,heatshock);

else error 'please enter plot_option 1 - 7'

end

