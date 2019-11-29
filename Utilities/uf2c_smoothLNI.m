function VetSmo = uf2c_smoothLNI(Vet,n)
%
%Linear Smooth Function
%
%VetSmo = smoothLNI(Vet,n)
%Vet = vector that you want to smooth
%n = size of the moving average
% 
%
%Brunno Machado de Campos 09/10/2012
%University of Campinas

%To contact me, highlight the next line (without the %) and press F9:
%email = [char(98) char(114) 'unnocampos' num2str(cos(2*pi)) char(64) char(116) char(101) 'rra.com.br']

Comp = length(Vet);

if mod(n,2) == 0
    n = n-1;
end
 
VetSmo = zeros(Comp,1);

for i = 1:Comp
    
    if i<=((n-1)/2)
        VetSmo(i) = mean(Vet(i-(i-1):i+(i-1)));
    else if i >= (Comp - ((n-1)/2))
        VetSmo(i) = mean(Vet(i-(Comp-i):i+(Comp-i)));    
        else    
    VetSmo(i) = mean(Vet(i-((n-1)/2):i+((n-1)/2)));    
         end
    end
end
