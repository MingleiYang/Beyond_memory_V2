function [counts] = lifetimefunction(state, W)
%lifetimefunction  calculates consecutive counts from state vector.
% Leonie Ringrose, 17.05.19

% Used for PRE/TRE - promoter model to calculate counts and transitions from state vector, this is used to calculate lifetimes. 


% States : PRE/TRE (1=M, -1=A, 0=U); PROMOTER (1=B, -1=F, 0= neither).
%---------------------------------------------------------------------------------
%
% counts column 1 = state 1 counts
% counts column 2 = state 1 transitions
% counts column 3 = state -1 counts
% counts column 4 = state -1 transitions
% counts column 5 = state 0 counts
% counts column 6 = state 0 transitions


counts = zeros (W,6); %set up output matrix.

% set up the first row of counts so that it reflects the correct state of
% the PRE/TRE or the promoter.
%-------------------------------------------------------------------------------

if state (1) == 1
    counts (1,1)= 1;
elseif state (1)== -1;
    counts (1,3)= 1;
else counts (1,5)=1;
end

% calculate counts 
%-------------------------------------------------------------------------------

for h = 1:W-1; 
    
   if state (h)==1; 
       
        if state (h+1) == 1
            counts (h+1,1)= counts (h,1)+1;
            counts (h+1,2)= counts (h,2);
            counts (h+1,3)= counts (h,3);
            counts (h+1,4)= counts (h,4);
            counts (h+1,5)= counts (h,5);
            counts (h+1,6)= counts (h,6);
        
        else
            counts (h+1,1)= counts (h,1);   % counts for state 1 no change 
            counts (h+1,2)= counts (h,2)+1; %transition away from state 1.
            
            if state (h+1) == -1
                counts (h+1,3)= counts (h,3)+1; % gain state -1
                counts (h+1,5)= counts (h,5); % no change state 0
            else 
                counts (h+1,5)= counts (h,5)+1; % gain state 0
                counts (h+1,3)= counts (h,3); % no change state -1
            end
                
            counts (h+1,4)= counts (h,4); % no transition away from state -1
            counts (h+1,6)= counts (h,6); % no transition away from state 0
        end 
        
   elseif state (h)==-1; 
       
        if state (h+1) == -1
            counts (h+1,1)= counts (h,1);
            counts (h+1,2)= counts (h,2);
            counts (h+1,3)= counts (h,3)+1;
            counts (h+1,4)= counts (h,4);
            counts (h+1,5)= counts (h,5);
            counts (h+1,6)= counts (h,6);
        
        else
            counts (h+1,3)= counts (h,3);   % counts for state -1 no change 
            counts (h+1,4)= counts (h,4)+1; %transition away from state -1.
            
            if state (h+1) == 1
                counts (h+1,1)= counts (h,1)+1; % gain state 1
                counts (h+1,5)= counts (h,5); % no change state 0
            else 
                counts (h+1,5)= counts (h,5)+1; % gain state 0
                counts (h+1,1)= counts (h,1); % no change state 1
            end
                
            counts (h+1,2)= counts (h,2); % no transition away from state 1
            counts (h+1,6)= counts (h,6); % no transition away from state 0
        end 
       
   else % state (h) ==0
        if state (h+1) == 0
            counts (h+1,1)= counts (h,1);
            counts (h+1,2)= counts (h,2);
            counts (h+1,3)= counts (h,3);
            counts (h+1,4)= counts (h,4);
            counts (h+1,5)= counts (h,5)+1;
            counts (h+1,6)= counts (h,6);
        
        else
            counts (h+1,5)= counts (h,5);   % counts for state 0 no change 
            counts (h+1,6)= counts (h,6)+1; %transition away from state 0.
            
            if state (h+1) == 1
                counts (h+1,1)= counts (h,1)+1; % gain state 1
                counts (h+1,3)= counts (h,3); % no change state -1
            else 
                counts (h+1,3)= counts (h,3)+1; % gain state -1
                counts (h+1,1)= counts (h,1); % no change state 1
            end
                
            counts (h+1,2)= counts (h,2); % no transition away from state 1
            counts (h+1,4)= counts (h,4); % no transition away from state -1
        end 
   end
   

    
end % h loop (counts in state matrix)

counts (W,:) = counts (W-1,:); % takes the last row and copies it.



 
end

