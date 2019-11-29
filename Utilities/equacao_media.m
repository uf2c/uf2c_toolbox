%gera valores para a equação de template ou médias em geral no incalc
%colocar o numero de imagens   'n'
%o resultado será impresso na command window

n=209;


b='i1';
for j=2:n
a=sprintf('%s+i%d',b,j);
b=a;
end
sprintf('(%s)/%d',b,j)