 function [p3input,p4input,p3_3input,p4_3input,minput,p5input,Cp_einput1,Ce_pinput1,Cp_einput2,Ce_pinput2,Cp_einput3,Ce_pinput3,tfinput,parameter2input] = memory_inputs(parameter1,parameter2,p3ser,p4ser,p3_3ser,p4_3ser,mser,p5ser,Cp_eser1,Ce_pser1,Cp_eser2,Ce_pser2,Cp_eser3,Ce_pser3,tfser);

% memory_inputs   calculates input vectors for P3Val et al. 
% Leonie Ringrose, 17.05.19

% used for parameter space plots.

x= numel(parameter1); 
y= numel (parameter2);

%parameter 1. This is on the y axis of parameter space plots 

p3input= p3ser;
p4input= p4ser;
p3_3input= p3_3ser;
p4_3input= p4_3ser;
p5input= p5ser;
minput= mser;
Cp_einput1 = Cp_eser1;
Ce_pinput1 = Ce_pser1;
Cp_einput2 = Cp_eser2;
Ce_pinput2 = Ce_pser2;
Cp_einput3 = Cp_eser3;
Ce_pinput3  = Ce_pser3;
tfinput = tfser;

if y > 1
    
    for i = 1:y-1 % Adds each series  to itself y times.

        p3input = [p3input,p3ser]; 
        p4input = [p4input,p4ser];
        p3_3input = [p3_3input,p3_3ser];  
        p4_3input = [p4_3input,p4_3ser];
        p5input = [p5input,p5ser];
        minput = [minput,mser];
        Cp_einput1 = [Cp_einput1,Cp_eser1];
        Ce_pinput1 = [Ce_pinput1,Ce_pser1];
        Cp_einput2 = [Cp_einput2,Cp_eser2];
        Ce_pinput2 = [Ce_pinput2,Ce_pser2];
        Cp_einput3 = [Cp_einput3,Cp_eser3];
        Ce_pinput3  = [Ce_pinput3,Ce_pser3];
        tfinput = [tfinput,tfser];

    end

end

%parameter 2. Now redefine the vector for parameter 2. This is on the x axis of parameter space plots 

  
    for i = 1:x % define top and bottom in terms of y, repeat x times.
    
    parameter2input(i*y-y+1:i*y) = parameter2(i);

    end


if parameter2 == p3ser
    p3input = parameter2input;

end

if parameter2 == p3_3ser
    p3_3input = parameter2input;

end

if parameter2 == p4ser
    p4input = parameter2input;

end
if parameter2 == p4_3ser
    p4_3input = parameter2input;

end

if parameter2 == p5ser
    p5input = parameter2input;

end

if parameter2 == mser
    minput = parameter2input;

end

if parameter2 == Cp_eser1
    Cp_einput1 = parameter2input;
    Ce_pinput1 = parameter2input;
end


if parameter2 == Cp_eser2
    Cp_einput2 = parameter2input;
    Ce_pinput2 = parameter2input;

end


if parameter2 == Cp_eser3
    Cp_einput3 = parameter2input;
    Ce_pinput3 = parameter2input;
end


if parameter2 == tfser
    tfinput = parameter2input;

end


clear x
clear y
clear i

end

