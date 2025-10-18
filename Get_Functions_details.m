% The 23 test functions are from:
% Xin Yao, Yong Liu, & Guangming Lin. (1999). Evolutionary programming made faster. 
% IEEE Transac-tions on Evolutionary Computation, 3(2), 82–102. https://doi.org/10.1109/4235.771163

% The three-bar truss design problem is from:
% RAY, T., & SAINI, P. (2001). Engineering design optimization using a swarm with an intelligent infor-mation sharing among individuals. 
% Engineering Optimization, 33(6), 735–748. https://doi.org/10.1080/03052150108940941

% The pressure vessel design problem is from:
% Wang, Z., Luo, Q., & Zhou, Y. (2020). Hybrid metaheuristic algorithm using butterfly and flower polli-nation base on mutualism mechanism for global optimization problems. 
% Engineering with Computers, 37(4), 3665–3698. https://doi.org/10.1007/s00366-020-01025-8

% The speed reducer design problem is from:
% Zhou, Y., Zhang, S., Luo, Q., & Abdel-Baset, M. (2019). CCEO: Cultural cognitive evolution optimization algorithm. 
% Soft Computing, 23(23), 12561–12583. https://doi.org/10.1007/s00500-019-03806-w 

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=30;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=30;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=30;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=30;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
%         dim=30;
        dim=30;
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=30;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F13'
        fobj = @F13;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F14'
        fobj = @F14;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F15'
        fobj = @F15;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F16'
        fobj = @F16;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F17'
        fobj = @F17;
        lb=[-5,0];
        ub=[10,15];
        dim=2;
        
    case 'F18'
        fobj = @F18;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F19'
        fobj = @F19;
        lb=0;
        ub=1;
        dim=3;
        
    case 'F20'
        fobj = @F20;
        lb=0;
        ub=1;
        dim=6;     
        
    case 'F21'
        fobj = @F21;
        lb=0;
        ub=10;
        dim=4;    
%         dim=4;
    case 'F22'
        fobj = @F22;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F23'
        fobj = @F23;
        lb=0;
        ub=10;
        dim=4;
    
    case 'three_bar'
        fobj = @three_bar;
        lb=0;
        ub=1;
        dim=2;

    case 'pressure_vessel'
        fobj = @pressure_vessel;
        lb=[0,0,10,10];
        ub=[99,99,200,200];
        dim=4;
    case 'reducer_object'
        fobj = @reducer_object;
        lb=[2.6,0.7,17,7.3,7.3,2.9,5.0];
        ub=[3.6,0.8,28,8.3,8.3,3.9,5.5];
        dim=7;

end

end

% F1

function o = F1(x)
o=sum(x.^2);
end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4

function o = F4(x)
o=max(abs(x));
end

% F5

function o = F5(x)
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
dim=size(x,2);
o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];

for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F20

function o = F20(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F21

function o = F21(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F22

function o = F22(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F23

function o = F23(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

function o = three_bar(x)
x1=x(1); x2=x(2);
l=100; P=2; sigma=2;
f = (2 * sqrt(2) * x1 + x2) * l;
% Calculation of constraints
g = [
    (sqrt(2)*x1 + x2) / (sqrt(2)*x1^2 + 2*x1*x2) * P - sigma;
    (x2) / (sqrt(2)*x1^2 + 2*x1*x2) * P - sigma;
    (1) / (sqrt(2)*x2 + x1) * P - sigma;
];

violations = max(0, g);  % All constraint violations
penalty = sum(40 * exp(10.0 * violations) - 40);
o = f + penalty;
end

function o = pressure_vessel(x)
x1=x(1); x2=x(2); x3=x(3); x4=x(4);
f = 0.6224*x1*x3*x4 + 1.7781*x2*x3^2 + 3.1661*x1^2*x4 + 19.84*x1^2*x3;
% Calculation of constraints
g = [
    -x1 + 0.0193*x3;
    -x2 + 0.00954*x3;
    -pi*x3^2*x4 - (4/3)*pi*x3^3 + 1296000;
    x4 - 240;
];

violations = max(0, g);  % All constraint violations
penalty = sum(800 * exp(10.0 * violations) - 800);
o = f + penalty;
end

function o = reducer_object(x)
z1 = x(1); z2 = x(2); z3 = x(3); z4 = x(4);
z5 = x(5); z6 = x(6); z7 = x(7);
% Calculate the value of the objective function
f = 0.7894 * z1 * z2^2 * (3.3333 * z3^2 + 14.9334 * z3 - 43.0934) ...
    - 1.508 * z1 * (z6^2 + z7^2) ...
    + 0.7854 * (z5 * z7^2 + z4 * z6^2) ...
    + 7.477 * (z6^3 + z7^3);
% Calculation of constraints
g = [
    -z1 * z2^2 * z3 + 27;               % g1(x) ≤ 0 
    -z1 * z2^2 * z3^2 + 397.5;           % g2(x) ≤ 0
    -z2 * z6^4 * z3 * z4^-3 + 1.93;      % g3(x) ≤ 0
    -z2 * z7^4 * z3 * z5^-3 + 1.93;      % g4(x) ≤ 0
    10 * z6^-3 * sqrt((745*z4/(z2*z3))^2 + 16.9e6) - 1100;  % g5(x) ≤ 0
    10 * z7^-3 * sqrt((745*z5/(z2*z3))^2 + 157.5e6) - 850;  % g6(x) ≤ 0
    z2 * z3 - 40;                        % g7(x) ≤ 0
    -z1/z2 + 5;                          % g8(x) ≤ 0
    z1/z2 - 12;                          % g9(x) ≤ 0
    1.5*z6 - z4 + 1.9;                   % g10(x) ≤ 0
    1.1*z7 - z5 + 1.9                    % g11(x) ≤ 0
];

violations = max(0, g);  % All constraint violations
penalty = sum(40 * exp(10.0 * violations) - 40);
o = f + penalty;
end


function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));



end