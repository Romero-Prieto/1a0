function n = rOMaN(n)

if n < 4000
    Th = floor(n/1000);
    Hu = floor((n - Th*1000)/100);
    Te = floor((n - Th*1000 - Hu*100)/10);
    Di = (n - Th*1000 - Hu*100 - Te*10);
    
    Th = string(char('M'*ones(1,Th)));
    
    if mod(Hu,5) <= 3
        Hu = string(char('D'*ones(1,floor(Hu/5)))) + string(char('C'*ones(1,mod(Hu,5))));
    else
        Hu = "C" + string(char('M'*ones(1,floor(Hu/5)))) + string(char('D'*ones(1,ceil((-Hu + 5)/5))));
    end
    
    if mod(Te,5) <= 3
        Te = string(char('L'*ones(1,floor(Te/5)))) + string(char('X'*ones(1,mod(Te,5))));
    else
        Te = "X" + string(char('C'*ones(1,floor(Te/5)))) + string(char('L'*ones(1,ceil((-Te + 5)/5))));
    end
    
    if mod(Di,5) <= 3
        Di = string(char('V'*ones(1,floor(Di/5)))) + string(char('I'*ones(1,mod(Di,5))));
    else
        Di = "I" + string(char('X'*ones(1,floor(Di/5)))) + string(char('V'*ones(1,ceil((-Di + 5)/5))));
    end
    
    n  = char(Th + Hu + Te + Di);
end
