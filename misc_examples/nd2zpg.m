function  [zero pole gain ] = nd2zpg(num,den)
    den2 = den./den(1); 
    num2 = num./den(1); 
    gain = num2(1);
    zero = roots(num2)'; 
    pole = roots(den2)'; 
end  

