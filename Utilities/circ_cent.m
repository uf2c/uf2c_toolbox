function [C]=circ_cent(P1,P2,R)

%Function used to compute the center of the circle given two points and a
%radius. Takes only real inputs and gives only real outputs. It can
%modified easily to accomodate complex numbers as well. Remove the if
%condition for checking if inputs are real numbers and also the if
%condition to check if center has real value.
%Make sure the units are consistent

%Inputs:
    %P1-First point e.g. [3 4] => (x1,y1)
    %P2-Second point  
    %R-Radius
%Output:
    %C-Center of the circle. Its a 2x2 matrix. First row represents first 
    %possible center (x1,y1) and second row is the second possible center
    
    if(~(isreal(P1) && isreal(P2) && isreal(R)))
        error('Invalid/Complex inputs');
    end
    
        

a=P1(1)^2-P2(1)^2;
b=P1(2)^2-P2(2)^2;
c=-2*(P1(1)-P2(1));
d=-2*(P1(2)-P2(2));
e=a+b;
Coeff1=1+(d/c)^2;
Coeff2=((2*d*e/c^2)+(P2(1)*2*d/c)-2*P2(2));
Coeff3=((e/c)^2+(2*P2(1)*e/c)+P2(1)^2+P2(2)^2-R^2);
All_coeff=[Coeff1 Coeff2 Coeff3];
Eq_root=roots(All_coeff);
C=zeros(length(Eq_root),2);
for i=1:length(Eq_root)
    x=-(e+d*Eq_root(i))/c;
    C(i,1)=x;
    C(i,2)=Eq_root(i);
end

if(~isreal(C))
    error('Circle with the specified radius does not pass through specified points');
end
    