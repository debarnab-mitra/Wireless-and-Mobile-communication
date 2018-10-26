
clear;


//system specifications
n_points = 10000;

A = linspace(0.1,10,20)*(10^-3);

j = 1;

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = sqrt(10^(-PLdB/10)/2);



//outer loop for different amplitude values
while j <=20;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    

X = A(j)*complex(I,Q);

//fading
HR = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI = grand(n_points, 1, "nor", 0, 1/sqrt(2));


h = complex(HR,HI);

fc = 800*10^6;
i = 1;


//Noise
noise_powerdBm = -100;
sigma = sqrt(10^(-13)/2);
NR = grand(n_points, 1, "nor", 0, sigma);
NI = grand(n_points, 1, "nor", 0, sigma);

N = complex(NR,NI);
num_error = 0;


//inner loop for sending 10000 QPSK symbols
while i <= n_points

//signal power calculation
signal_powerdB = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI(i)^2+HR(i)^2);

if signal_powerdB < -165 then
    num_error = num_error + 1;
end

i = i+1;

end


Pe(j) = num_error/n_points;


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'bs-');
xtitle( 'Poutage vs Tx SNR', 'Tx SNR', 'Poutage')

