clear

//system specifications
f_samp = 2000*10^6;
fc = 800*10^6;
Ts = 1/f_samp;
n = linspace(1,2000,2000);


n_points = 10000;

A = linspace(0.1,10,20)*(10^-3);

j = 1;

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = sqrt(10^(-PLdB/10)/2);



//outer loop for different A values
while j <=20;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    
//coef of cos2pifct and sin2pifct
X = A(j)*complex(I,Q);

//fading
//3 independent fading values

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
No = 10^-13;
sigma = sqrt(No/2);
NR = grand(n_points, 1, "nor", 0, sigma);
NI = grand(n_points, 1, "nor", 0, sigma);

N = complex(NR,NI);
num_error = 0;

while i <= n_points

//y = |h|^2X + h*N

signal_powerdB_1 = HI_1(i)^2+HR_1(i)^2;
signal_powerdB_2 = (HI_2(i)^2+HR_2(i)^2);
signal_powerdB_3 = (HI_3(i)^2+HR_3(i)^2);

[index,val] = max(signal_powerdB_1,signal_powerdB_2,signal_powerdB_3);

//modulation step where only the channel with max SNR is taken
y = val*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);

//demodulation steps


    
    y_demod_I =  y*(cos(2*%pi*fc*n*Ts))';
    y_demod_Q =  y*(sin(2*%pi*fc*n*Ts))';
    




I_demod(i) = sign(y_demod_I);
Q_demod(i) = sign(y_demod_Q);

if (I(i) ~= I_demod(i)) or (Q(i) ~= Q_demod(i))
    num_error = num_error + 1;    
end

i = i+1;

end


Pe(j) = num_error/n_points;


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe));
xtitle( 'Pe vs Tx SNR, antenna diversity', 'Tx SNR', 'Pe')
