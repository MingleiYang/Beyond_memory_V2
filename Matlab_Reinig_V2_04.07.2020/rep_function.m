function [nucstart, sitestart] = rep_function(m,nucs,sites,ti,wipe)

% rep_function   Replication after timecourse.
% Leonie Ringrose, 17.05.19
% replicates the last row of nucs and sites and passes it back as the first row to timecourse function.


R= rand ([m,1]);%makes a list of m random numbers (one for each nucleosome)
    Rnucs= nucs (ti,:);%extracts the last row of nucs 
    Rsites = sites (ti,:);% extracts the last row of sites. In the current version all sites are set to 0 after rep.

    for r = 1:m  %one for each nucleosome (1-40)
        if  R(r)<=wipe % checks each random number in R at position 1-35. Should be 0.5!!
            Rnucs (r)=0; % sets the nucleosome at same position to zero in Rnucs.
        
        else
            Rnucs (r)= nucs (ti,r); %if not, it leaves it as it was.
        
        end

    end

    Rsites = 0; 
    
    nucstart = Rnucs;
    sitestart = Rsites;
    
end

