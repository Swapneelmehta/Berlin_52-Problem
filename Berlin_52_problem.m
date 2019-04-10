clc;
clear;
close all;

Berlin =[1         565         575
    2          25         185
    3         345         750
    4         945         685
    5         845         655
    6         880         660
    7          25         230
    8         525        1000
    9         580        1175
    10         650        1130
    11        1605         620
    12        1220         580
    13        1465         200
    14        1530           5
    15         845         680
    16         725         370
    17         145         665
    18         415         635
    19         510         875
    20         560         365
    21         300         465
    22         520         585
    23         480         415
    24         835         625
    25         975         580
    26        1215         245
    27        1320         315
    28        1250         400
    29         660         180
    30         410         250
    31         420         555
    32         575         665
    33        1150        1160
    34         700         580
    35         685         595
    36         685         610
    37         770         610
    38         795         645
    39         720         635
    40         760         650
    41         475         960
    42          95         260
    43         875         920
    44         700         500
    45         555         815
    46         830         485
    47        1170          65
    48         830         610
    49         605         625
    50         595         360
    51        1340         725
    52        1740         245 ];
x = Berlin(:,2);
y = Berlin(:,3);
n=prod(size(x));
Max_Iter = 30; % Maximum Number of Iterations
Na = 50;  
Dist = zeros(n,n);
Q = 100;
Tau_0 = 10^-6;
alpha = 1;              
beta = 5;
rho = 0.5;
T = Tau_0*ones(n,n);
anttour = zeros(Max_Iter,n);
antcost = zeros(n,1);
Bestcost = zeros(Max_Iter,1);
Bestcostsoln = inf;
for i = 1:n
    for j = 1:n
        
        Dist(i,j) = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
        
    end
end
eta = 1./Dist;
for iter = 1:Max_Iter
    for l = 1:n-1
    for k = 1:Na 
    
        anttour(k,1) = 1;
        i = anttour(k,l);
        P = (T(i,:).^alpha)./(Dist(i,:).^beta);
        P(anttour(k,1:l))=0;
        P = P/sum(P);
        r = rand;
        C = cumsum(P);
        j = find(r<=C,1,'first');
        anttour(k,l+1)=j;
    end
    end
    for k = 1:Na
        tour = [anttour(k,:) anttour(k,1)];
        L = 0;
        for i = 1:length(tour)-1
             L = L + Dist(tour(i),tour(i+1));
        end
        antcost(k,1) = L;
        if antcost(k,1) < Bestcostsoln
            Bestcostsoln = antcost(k,1);
            Besttoursoln = anttour(k,:);
        end
    end
    
    T = (1-rho)*T;
    for k = 1:Na
        tour = [anttour(k,:) anttour(k,1)];
        for l = 1:n
            i = tour(l);
            j = tour(l+1);
            T(i,j) = T(i,j)+Q/antcost(k,1);
        end
    end
    Bestcost(iter) = Bestcostsoln;
end
x_1 = ones(53,1);
y_1 = ones(53,1);
for i = 1:52
    x_1(i,1) = Berlin(Besttoursoln(1,i),2);
end
x_1(53,1) = x_1(1,1);
for i = 1:52
    y_1(i,1) = Berlin(Besttoursoln(1,i),3);
end
y_1(53,1) = y_1(1,1);
plot(x_1,y_1,'g-o');
figure
plot(1:Max_Iter,Bestcost)