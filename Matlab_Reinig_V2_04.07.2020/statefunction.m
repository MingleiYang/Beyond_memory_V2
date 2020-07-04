function [WA, WM, WB, WF, stateP, stateE] = statefunction (windowsize,t,step,W,Arun,Mrun,Brun,Frun)

%statefunction   calculates PRE/TRE and promoter states.

% Leonie Ringrose, 17.05.19

wpositions = (1:step:t);

WA = zeros (W,1);
WM = zeros (W,1);
WB = zeros (W,1);
WF = zeros (W,1);

stateP = zeros (W,1);
stateE = zeros (W,1);

for g = 1:W % sliding window for states.
    
    pos = wpositions (g);
    
    WA(g,1) = mean(Arun(pos:pos+windowsize-1));
    WM(g,1) = mean(Mrun(pos:pos+windowsize-1));
    WB(g,1) = mean(Brun(pos:pos+windowsize-1));
    WF(g,1) = mean(Frun(pos:pos+windowsize-1));
    
    
    if WM(g,1) > 1.5* WA(g,1)
        stateP(g,1)= -1;          %assigns the state M if M >1.5 * A
    elseif WA(g,1) > 1.5*WM(g,1)  %assigns the state A if A >1.5 * M
        stateP(g,1) = 1;
    else stateP(g,1)= 0;          %assigns the state U if neither is true
    end
    
    if WB(g,1) > 1.5* WF(g,1)     %assigns the state B if B >1.5 * F
        stateE(g,1)= 1;
    elseif WF(g,1) > 1.5*WB(g,1)  %assigns the state F if F >1.5 * B
        stateE(g,1) = -1;
    else stateE(g,1)= 0;          %assigns the state N if neither is true
    end
   
end % g loop (sliding window for states)

end

