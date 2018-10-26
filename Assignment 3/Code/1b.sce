clear;

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



// outer loop for different amplitude values
while j <=20;
    
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;    
//coef of cos2pifct and sin2pifct
X = A(j)*complex(I,Q);

//fading
HR = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI = grand(n_points, 1, "nor", 0, 1/sqrt(2));

H_sq = grand(n_points,1,"exp",1/sqrt(2));

h = complex(HR,HI);

fc = 800*10^6;
i = 1;


//Noise
noise_powerdBm = -100;
//sigma = (0.001*(10^(noise_powerdBm/20)))/2;

No = 10^-13;
sigma = sqrt(No/2);
NR = grand(n_points, 1, "nor", 0, sigma);
NI = grand(n_points, 1, "nor", 0, sigma);

N = complex(NR,NI);
num_error = 0;
bit_error = 0;

//inner loop for sending 10000 QPSK symbols
while i <= n_points

//for each X take 2000 samples of cos2pifct and sine2pifct and multiply with X
//then multiple with PL*h
//then addd noise
//then again multiple with same samples and avg

//this is the modulation step
    y = (HR(i)^2 + HI(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);
    



//demodulation steps // corelator receiver kind

    y_demod_I =  y*(cos(2*%pi*fc*n*Ts))';
    y_demod_Q =  y*(sin(2*%pi*fc*n*Ts))';
    




I_demod(i) = sign(y_demod_I);
Q_demod(i) = sign(y_demod_Q);


//symbol error calculation
if (I(i) ~= I_demod(i)) or (Q(i) ~= Q_demod(i))
    num_error = num_error + 1;    
end

//if (I(i) ~= I_demod(i))
//    bit_error = bit_error + 1;    
//end

//if (Q(i) ~= Q_demod(i))
 //   bit_error = bit_error + 1;    
//end

i = i+1;

end


Pe(j) = num_error/n_points;
P_bit_error(j) = bit_error/(2*n_points);

j = j+1;
end



//sym_t = 10^-6;

//BER = P_bit_error/sym_t;

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'bs-');
//plot(10*log10((A.^2)/(10^-13)),(P_bit_error),'bs-');
//plot(10*log10((A.^2)/(10^-13)),BER,'bs-');
xtitle( 'Pe vs Tx SNR', 'Tx SNR', 'Pe')
