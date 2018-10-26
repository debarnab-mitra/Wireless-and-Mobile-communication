

clear;

n_points = 10000;

A = linspace(0.1,10,40)*(10^-3);

j = 1;

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = sqrt(10^(-PLdB/10)/2);




while j <=40;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    
//coef of cos2pifct and sin2pifct
X = A(j)*complex(I,Q);

//fading
HR = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI = grand(n_points, 1, "nor", 0, 1/sqrt(2));

//H_sq = grand(n_points,1,"exp",1/sqrt(2));

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

while i <= n_points

//y = |h|^2X + h*N
//y = (HR(i)^2 + HI(i)^2)*PL*X(i) + N(i);
//y = (sqrt(H_sq(i)))*PL*X(i) + N(i);


//I_demod(i) = sign(real(y));
//Q_demod(i) = sign(imag(y));

signal_powerdB = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI(i)^2+HR(i)^2);

if signal_powerdB < -165 then
    num_error = num_error + 1;
end

i = i+1;

end


Pe(j) = num_error/n_points;


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','blue');
//xtitle( 'Poutage vs Tx SNR', 'Tx SNR', 'Poutage')

//////////////////////////////////////////////////////////////////////////////////////////////
clear;



n_points = 10000;

A = linspace(0.1,10,40)*(10^-3);

j = 1;

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = sqrt(10^(-PLdB/10)/2);




while j <=40;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    
//coef of cos2pifct and sin2pifct
X = A(j)*complex(I,Q);

//fading
HR_1 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_1 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
//H_sq = grand(n_points,1,"exp",1/sqrt(2));


HR_2 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_2 = grand(n_points, 1, "nor", 0, 1/sqrt(2));

HR_3 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_3 = grand(n_points, 1, "nor", 0, 1/sqrt(2));

//H_sq = grand(n_points,1,"exp",1/sqrt(2));

h = complex(HR,HI);

fc = 800*10^6;
i = 1;


//Noise
noise_powerdBm = -100;
No = 10^-13;
sigma = sqrt(No/2);
NR = grand(n_points, 1, "nor", 0, sigma);
NI = grand(n_points, 1, "nor", 0, sigma);

N = complex(NR,NI);
num_error = 0;

while i <= n_points

//y = |h|^2X + h*N
//y = (HR(i)^2 + HI(i)^2)*PL*X(i) + N(i);
//y = (sqrt(H_sq(i)))*PL*X(i) + N(i);


//I_demod(i) = sign(real(y));
//Q_demod(i) = sign(imag(y));

signal_powerdB_1 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_1(i)^2+HR_1(i)^2);
signal_powerdB_2 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_2(i)^2+HR_2(i)^2);
signal_powerdB_3 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_3(i)^2+HR_3(i)^2);

signal_powerdB = max(signal_powerdB_1,signal_powerdB_2,signal_powerdB_3);

if signal_powerdB < -165 then
    num_error = num_error + 1;
end

i = i+1;

end


Pe(j) = num_error/n_points;


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','red');
//xtitle( 'Poutage vs Tx SNR , antenna diversity', 'Tx SNR', 'Poutage')
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

clear;

n_points = 10000;

A = linspace(0.1,10,40)*(10^-3);

j = 1;

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = sqrt(10^(-PLdB/10)/2);




while j <=40;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    
//coef of cos2pifct and sin2pifct
X = A(j)*complex(I,Q);

//fading
HR = zeros(3,n_points);
HI = zeros(3,n_points);

HR_1 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_1 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
//H_sq = grand(n_points,1,"exp",1/sqrt(2));


HR_2 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_2 = grand(n_points, 1, "nor", 0, 1/sqrt(2));

HR_3 = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI_3 = grand(n_points, 1, "nor", 0, 1/sqrt(2));



h = complex(HR,HI);

fc = 800*10^6;
i = 1;


//Noise
noise_powerdBm = -100;
sigma = (0.001*(10^(noise_powerdBm/20)))/2;
NR_1 = grand(n_points, 1, "nor", 0, sigma);
NI_1 = grand(n_points, 1, "nor", 0, sigma);

NR_2 = grand(n_points, 1, "nor", 0, sigma);
NI_2 = grand(n_points, 1, "nor", 0, sigma);

NR_3 = grand(n_points, 1, "nor", 0, sigma);
NI_3 = grand(n_points, 1, "nor", 0, sigma);

//N1 = complex(NR,NI);
//N2 = complex(NR,NI);
//N3 = complex(NR,NI);
num_error = 0;

while i <= n_points

//y = |h|^2X + h*N
//y = (HR(i)^2 + HI(i)^2)*PL*X(i) + N(i);
//y = (sqrt(H_sq(i)))*PL*X(i) + N(i);


//I_demod(i) = sign(real(y));
//Q_demod(i) = sign(imag(y));



signal_powerdB_1 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_1(i)^2+HR_1(i)^2);
signal_powerdB_2 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_2(i)^2+HR_2(i)^2);
signal_powerdB_3 = 10*log10(2*A(j)^2) - PLdB + 10*log10(HI_3(i)^2+HR_3(i)^2);

//signal_powerdB = max(signal_powerdB_1,signal_powerdB_2,signal_powerdB_3);

a = 0;

if signal_powerdB_1 < -165 then
    a = a + 1;
end

if signal_powerdB_2 < -165 then
    a = a + 1;
end

if signal_powerdB_3 < -165 then
    a = a + 1;
end

if(a >= 2)
   
   num_error = num_error + a; 
end    


i = i+1;

end


Pe(j) = num_error/(3*n_points);


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','green');
xtitle( 'Poutage vs Tx SNR', 'Tx SNR', 'Poutage');
legend(["No diversity","Antenna diversity","Time diversity"]);

