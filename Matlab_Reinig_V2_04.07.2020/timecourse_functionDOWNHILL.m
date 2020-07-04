function [nucs, sites, OUT,DOWNHILLOUT] = timecourse_functionDOWNHILL (C_down,h_down,ti,m,tf,model,p3,p4,p5,p6,Cp_e,Ce_p,p1,p2,nucstart,sitestart,couple_ep_p3start,couple_ep_p4start)
  
% timecourse_functionDOWNHILL  simulates PRE/TRE and promoter over time t.
% Leonie Ringrose, 17.05.19
% Specialised to adjust p1 via Hill function for eya gradient decay in eye disc

% set up output matrices 
%----------------------------------------------------------------------------------------------------------
nucs = zeros (ti,m); 
sites = zeros (ti,tf); 
OUT = zeros (ti,18); 
DOWNHILLOUT = zeros (ti,1);

% set up first row
%----------------------------------------------------------------------------------------------------------
 nucs (1,:)= nucstart; % sets all nucs in the first row to the same value. 
 sites(1,:)= sitestart; % sets all sites in the first row to the same value. Default = 0 
 
 couple_ep_p3= couple_ep_p3start; % set to 1 for first iteration, is then calculated from the enhancer state. 
 couple_ep_p4= couple_ep_p4start; %
    
    
for i= 1:ti
   

    %modify PRE
    %--------------------------------------------------------------------------------------------------------
    
    x = randi (m,1); % generates a random integer between 1 and the number of nucleosomes, this defines n1. 
    y = randi (m,1); % generates a second random integer, this defines n2.(if n1 and n2 are the same, nothing happens)

    r3= rand;
    r4= rand;
    r5= rand;
    r6= rand;
    
    
    X=nucs (i,:);   %extracts the ith row of nucs.  
    P=X(X>0);
    N=X(X<0);
    Z=X(X==0);      % Extracts +1, -1 and 0.
    M=numel(P);
    A=numel(N);
    U=numel(Z);     % counts them.
    
    OUT (i,1)= M;
    OUT (i,2)= A;
    OUT (i,3)= U; 
    OUT (i,4)= M/m; 
    OUT (i,5)= A/m; % puts them in OUT.
    
    nucs(i+1,:) = nucs(i,:);
    
    %calculate p3co, p4co
    %------------------------------------------------------------------------------------------------------------------------------------------
    
    p3co= couple_ep_p3*p3; %calculates p3co based on the value of couple_ep_p3. This is first set to 1, then calculated from enhancer output.
    p4co= couple_ep_p4*p4; %calculates p4co based on the value of couple_ep_p4. This is first set to 1, then calculated from enhancer output.
    
   
    %attempt recruited conversion to M
    %----------------------------------------------------------------------------------------------------------------------------------------------------
    
    if nucs(i,y)==1 % if n2 is in M state
       
        if r3<=p3co   % attempt conversion towards M state with probability p3co.        
            if nucs(i,y)> nucs(i,x)   % if n1 is not already in M state (ie it is 0 or -1)         
            nucs(i+1,x)= nucs(i,x)+1;  % modify n1 towards M. 
            else nucs(i+1,x)= nucs(i,x); % If n1 is already in M state, do nothing. 
            end      
        else nucs(i+1,x)= nucs(i,x); % if p3co not satisfied, do nothing.  
        end
    else nucs(i+1,x)= nucs(i,x); % if n2 not in M state, do nothing.
    
    end
                
        a = nucs(i+1,x); % store the modified or unmodified nucleosome as "a".
    
    %attempt recruited conversion to A
    %----------------------------------------------------------------------------------------------------
        
    if nucs(i,y)==-1 % if n2 is in A state
       
        if r4<=p4co   % attempt conversion towards A state with probability p4co.        
            if nucs(i,y)< nucs(i,x)   % if n1 is not already in A state (ie it is 1 or 0)         
            nucs(i+1,x)= nucs(i,x)-1;  % modify n1 towards A. 
            else nucs(i+1,x)= a; % If n1 is already in A state, leave it as previous value. 
            end      
        else nucs(i+1,x)= a; % if p4co not satified, do nothing.  
        end
    else nucs(i+1,x)= a; % if n2 not in A state, do nothing.
    
    end
    
    b = nucs(i+1,x); % store the modified or unmodified nucleosome as "b".
  
    %attempt noisy conversion on the modified or unmodified nucleosome.  
    %----------------------------------------------------------------------------------------------------
    if (r5<=p5) % probability of noisy conversion
          if b == 1 
            nucs(i+1,x) = b-1;  % if at 1, go down
            
          elseif b == -1
            nucs(i+1,x) = b+1;  % if at -1, go up
            
          else % if at 0, go up or down with probability p6. p6 is 0.5.
                
                if (r6<=p6)
                    nucs(i+1,x) = b+1;
                else
                    nucs(i+1,x) = b-1;
                end               
          end                            
    else nucs(i+1,x)= b;
    
    end
  
    %modify the enhancer. 
    %----------------------------------------------------------------------------------------------------
    
    z = randi (tf,1); % generates a random integer between 1 and the number of sites, this defines s1. 
    r1= rand;
    r2= rand;
   
    E=sites(i,:);   %extracts the ith row of nucs.  
    P=E(E==1);
    Z=E(E==0);      % Extracts +1 and 0.
    B=numel(P);     % Number of bound sites, B
    F=numel(Z);     % Number of Free sites, F.
    
    OUT (i,7)= B;   % puts them in OUT
    OUT (i,8)= F;   
    OUT (i,9)= B/tf;
    OUT (i,10)= F/tf;
    
    
    % Calculate the coupling parameters according to the model being used.
    % p1 and p2 affect the enhancer and depend on the PRE
    % p3 and p4 affect the PRE and depend on the enhancer. 
    %---------------------------------------------------------------------------------------------------------
   
    
    
     if model == 1 % p1, p2, p3, p4. Model 1 in paper. 
        couple_p1 = 2*(exp(Cp_e*(1/m)*(A-M)))/(1+(exp(Cp_e*(1/m)*(A-M))));%p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = 2*(exp(Cp_e*(1/m)*(M-A)))/(1+(exp(Cp_e*(1/m)*(M-A))));%p2 (enhancer transition to "F") scales with silent PRE.
        couple_ep_p3 = 2*(exp(Ce_p*(1/tf)*(F-B)))/(1+(exp(Ce_p*(1/tf)*(F-B)))); % p3 (PRE transition to "M") scales with silent enhancer.
        couple_ep_p4 = 2*(exp(Ce_p*(1/tf)*(B-F)))/(1+(exp(Ce_p*(1/tf)*(B-F))));% p4 (PRE transition to "A") scales with active enhancer.
    
    elseif model == 2 % p1, p2, p3, p4 Model 2 in paper. Needs C = 0.5* C for model 1.  
        couple_p1 = exp(Cp_e*(1/m)*(A-M));%p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = exp(Cp_e*(1/m)*(M-A));%p2 (enhancer transition to "F") scales with silent PRE.  
        couple_ep_p3 = exp(Ce_p*(1/tf)*(F-B)); % p3 (PRE transition to "M") scales with silent enhancer.     
        couple_ep_p4 = exp(Ce_p *(1/tf)*(B-F));% p4 (PRE transition to "A") scales with active enhancer.  
       
    elseif model == 3 % Model 1, p1, p4
        couple_p1 = 2*(exp(Cp_e*(1/m)*(A-M)))/(1+(exp(Cp_e*(1/m)*(A-M))));%p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = 1;
        couple_ep_p3 = 1;
        couple_ep_p4 = 2*(exp(Ce_p*(1/tf)*(B-F)))/(1+(exp(Ce_p*(1/tf)*(B-F))));% p4 (PRE transition to "A") scales with active enhancer.
    
    elseif model == 4 % Model 1, p2, p3
        couple_p1 = 1;
        couple_p2 = 2*(exp(Cp_e*(1/m)*(M-A)))/(1+(exp(Cp_e*(1/m)*(M-A))));%p2 (enhancer transition to "F") scales with silent PRE.
        couple_ep_p3 = 2*(exp(Ce_p*(1/tf)*(F-B)))/(1+(exp(Ce_p*(1/tf)*(F-B)))); % p3 (PRE transition to "M") scales with silent enhancer.
        couple_ep_p4 = 1;
        
    elseif model == 5 % Model 1, p1, p3
        couple_p1 = 2*(exp(Cp_e*(1/m)*(A-M)))/(1+(exp(Cp_e*(1/m)*(A-M))));%p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = 1;
        couple_ep_p3 = 2*(exp(Ce_p*(1/tf)*(F-B)))/(1+(exp(Ce_p*(1/tf)*(F-B)))); % p3 (PRE transition to "M") scales with silent enhancer.
        couple_ep_p4 = 1;
    
    elseif model == 6 % Model 1, p2, p4
        couple_p1 = 1;
        couple_p2 = 2*(exp(Cp_e*(1/m)*(M-A)))/(1+(exp(Cp_e*(1/m)*(M-A))));%p2 (enhancer transition to "F") scales with silent PRE.
        couple_ep_p3 = 1;
        couple_ep_p4 = 2*(exp(Ce_p*(1/tf)*(B-F)))/(1+(exp(Ce_p*(1/tf)*(B-F))));% p4 (PRE transition to "A") scales with active enhancer.
        
    elseif model ==  7 % Model 1, p1, p2, p3, p4 are only REDUCED. 
        couple_p1 = 2*(exp(Cp_e*(1/m)*(A-M)))/(1+(exp(Cp_e*(1/m)*(A-M))));%p1 (enhancer transition to "B") scales with active PRE. 
        
        if couple_p1 > 1
           couple_p1 = 1;
        end
 
        couple_p2 = 2*(exp(Cp_e*(1/m)*(M-A)))/(1+(exp(Cp_e*(1/m)*(M-A))));%p2 (enhancer transition to "F") scales with silent PRE.
        
        if couple_p2 > 1
           couple_p2 = 1;
        end
        
        couple_ep_p3 = 2*(exp(Ce_p*(1/tf)*(F-B)))/(1+(exp(Ce_p*(1/tf)*(F-B)))); % p3 (PRE transition to "M") scales with silent enhancer.
        
        if couple_ep_p3 > 1
           couple_ep_p3 = 1;
        end
        
        couple_ep_p4 = 2*(exp(Ce_p*(1/tf)*(B-F)))/(1+(exp(Ce_p*(1/tf)*(B-F))));% p4 (PRE transition to "A") scales with active enhancer.    
        
        if couple_ep_p4 > 1
           couple_ep_p4 = 1;
        end
        
    elseif model ==  8 % Model 1, p1, p2, p3, p4 are only INCREASED. 
        couple_p1 = 2*(exp(Cp_e*(1/m)*(A-M)))/(1+(exp(Cp_e*(1/m)*(A-M))));%p1 (enhancer transition to "B") scales with active PRE. 
        
            if  couple_p1 < 1
                couple_p1 = 1;
            end
 
        couple_p2 = 2*(exp(Cp_e*(1/m)*(M-A)))/(1+(exp(Cp_e*(1/m)*(M-A))));%p2 (enhancer transition to "F") scales with silent PRE.
        
            if  couple_p2 < 1
                couple_p2 = 1;
            end
        
        couple_ep_p3 = 2*(exp(Ce_p*(1/tf)*(F-B)))/(1+(exp(Ce_p*(1/tf)*(F-B)))); % p3 (PRE transition to "M") scales with silent enhancer.
        
            if  couple_ep_p3 < 1
                couple_ep_p3 = 1;
            end
        
        couple_ep_p4 = 2*(exp(Ce_p*(1/tf)*(B-F)))/(1+(exp(Ce_p*(1/tf)*(B-F))));% p4 (PRE transition to "A") scales with active enhancer.    
        
            if  couple_ep_p4 < 1
                couple_ep_p4 = 1;
            end
            
    elseif model == 9 % Model 2, p1, p4
        couple_p1 = exp(Cp_e*(1/m)*(A-M)); %p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = 1;
        couple_ep_p3 = 1;
        couple_ep_p4 = exp(Ce_p *(1/tf)*(B-F));% p4 (PRE transition to "A") scales with active enhancer.
    
    elseif model == 10 % Model 2, p2, p3
        couple_p1 = 1;
        couple_p2 = exp(Cp_e*(1/m)*(M-A)); %p2 (enhancer transition to "F") scales with silent PRE.
        couple_ep_p3 = exp(Ce_p*(1/tf)*(F-B));% p3 (PRE transition to "M") scales with silent enhancer.
        couple_ep_p4 = 1;
        
    elseif model == 11 % Model 2, p1, p3
        couple_p1 = exp(Cp_e*(1/m)*(A-M)); %p1 (enhancer transition to "B") scales with active PRE. 
        couple_p2 = 1;
        couple_ep_p3 = exp(Ce_p*(1/tf)*(F-B)); % p3 (PRE transition to "M") scales with silent enhancer.
        couple_ep_p4 = 1;
    
    elseif model == 12 % Model 2, p2, p4
        couple_p1 = 1;
        couple_p2 = exp(Cp_e*(1/m)*(M-A));%p2 (enhancer transition to "F") scales with silent PRE.
        couple_ep_p3 = 1;
        couple_ep_p4 = exp(Ce_p *(1/tf)*(B-F)); % p4 (PRE transition to "A") scales with active enhancer.
        
    elseif model ==  13 % Model 2, p1, p2, p3, p4 are only REDUCED. 
        couple_p1 = exp(Cp_e*(1/m)*(A-M)); %p1 (enhancer transition to "B") scales with active PRE. 
        
        if couple_p1 > 1
           couple_p1 = 1;
        end
        
        couple_p2 = exp(Cp_e*(1/m)*(M-A));%p2 (enhancer transition to "F") scales with silent PRE.
        if couple_p2 > 1
           couple_p2 = 1;
        end
        
        couple_ep_p3 = exp(Ce_p*(1/tf)*(F-B)); % p3 (PRE transition to "M") scales with silent enhancer.
        if couple_ep_p3 > 1
           couple_ep_p3 = 1;
        end
        
        couple_ep_p4 = exp(Ce_p *(1/tf)*(B-F));% p4 (PRE transition to "A") scales with active enhancer.    
        if couple_ep_p4 > 1
           couple_ep_p4 = 1;
        end
        
    elseif model ==  14 % Model 2, p1, p2, p3, p4 are only INCREASED. 
        couple_p1 = exp(Cp_e*(1/m)*(A-M)); %p1 (enhancer transition to "B") scales with active PRE. 
        
            if  couple_p1 < 1
                couple_p1 = 1;
            end
        
        couple_p2 = exp(Cp_e*(1/m)*(M-A));%p2 (enhancer transition to "F") scales with silent PRE.
        
            if  couple_p2 < 1
                couple_p2 = 1;
            end
        
        couple_ep_p3 = exp(Ce_p*(1/tf)*(F-B)); % p3 (PRE transition to "M") scales with silent enhancer.
            
            if  couple_ep_p3 < 1
                couple_ep_p3 = 1;
            end
        
        couple_ep_p4 = exp(Ce_p *(1/tf)*(B-F));% p4 (PRE transition to "A") scales with active enhancer.    
            if  couple_ep_p4 < 1
                couple_ep_p4 = 1;
            end        
            
    else error ('please choose model 1-14')
    end
    
    OUT (i,11)= couple_p1;
    OUT (i,12)= couple_p2;
    OUT (i,13)= couple_ep_p3; %initially set to 1, store in OUT.
    OUT (i,14)= couple_ep_p4; %initially set to 1, store in OUT.
    
   % Calculate the DOWNHILL function: 
   
    Hill = 1-((C_down*(i/100)^h_down/(1+(i/100)^h_down)));
    DOWNHILLOUT (i) = Hill;
    
   
     % Calculate the value of p1 using the coupling parameter for the effect of the PRE on the enhancer, depending on PRE state. 
    %---------------------------------------------------------------------------------------------------------
    
    p1co=p1*couple_p1*Hill; % use for the PRE affecting the enhancer. 
    p2co=p2*couple_p2; % use for the PRE affecting the enhancer. 
    
    OUT (i,15)= p1co; %store p1co - p4co- these are the modified p1- p4 values with coupling
    OUT (i,16)= p2co; 
    OUT (i,17)= p3co;
    OUT (i,18)= p4co;
    
    
    sites(i+1,:) = sites(i,:);% this takes the value of the last row and gives it to the next update (which will change it).
    
    
    % Attempt binding at site s1 NB I changed the program slightly on 3rd
    % Feb 2015, to treat either 1 or 0 . 
    %---------------------------------------------------------------------------------------------------------
    
if sites(i,z)==0
    
    if r1<=p1co %attempt binding at s1 with probability p1co.
       sites (i+1,z) = 1; % it goes to 1.
    else
       sites(i+1,z) = sites(i,z);% if 0 it stays 0.
    end
    
else % if it is 1, attempt unbinding
        if r2<=p2co % attempt unbinding at s1 with probability p2.      
        sites(i+1,z) = 0;  % if 0 it stays 0, if 1 it goes to 0        
        else  
        sites(i+1,z)= sites(i,z); % if r2 > p2, no change, it stays bound 
        end                
end        
        
end % i loop, time course. 1:ti.


end

