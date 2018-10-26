clear
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
while i <= n_points

////y = |h|^2X + h*N
//y = (HR(i)^2 + HI(i)^2)*PL*X(i) + N(i);


//instead of this
//for each X take 2000 samples of cos2pifct and sine2pifct and multiply with X
//then multiple with PL*h
//then addd noise
//then again multiple with same samples and avg


    
    y = (HR(i)^2 + HI(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);
    



//demodulation steps


    
    y_demod_I =  y*(cos(2*%pi*fc*n*Ts))';
    y_demod_Q =  y*(sin(2*%pi*fc*n*Ts))';
    




I_demod(i) = sign(y_demod_I);
Q_demod(i) = sign(y_demod_Q);

if (I(i) ~= I_demod(i)) or (Q(i) ~= Q_demod(i))
    num_error = num_error + 1;    
end

if (I(i) ~= I_demod(i))
    bit_error = bit_error + 1;    
end

if (Q(i) ~= Q_demod(i))
    bit_error = bit_error + 1;    
end

i = i+1;

end


Pe(j) = num_error/n_points;
P_bit_error(j) = bit_error/(2*n_points);

j = j+1;
end



sym_t = 10^-6;

BER = P_bit_error/sym_t;

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','blue');
//plot(10*log10((A.^2)/(10^-13)),(P_bit_error),'bs-');
//plot(10*log10((A.^2)/(10^-13)),BER,'bs-');



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
clear
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




while j <=20;
    
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
//sigma = (0.001*(10^(noise_powerdBm/20)))/2;

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

//y = val*PL*X(i) + N(i);
//y = (sqrt(H_sq(i)))*PL*X(i) + N(i);

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

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','red');

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
clear
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




while j <=20;
    
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
//sigma = (0.001*(10^(noise_powerdBm/20)))/2;


No = 10^-13;
sigma = sqrt(No/2);

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


//y_1 = (HR_1(i)^2 + HI_1(i)^2)*PL*X(i) + N1(i);
//y_2 = (HR_2(i)^2 + HI_2(i)^2)*PL*X(i) + N2(i);
//y_3 = (HR_3(i)^2 + HI_3(i)^2)*PL*X(i) + N3(i);
//y = (sqrt(H_sq(i)))*PL*X(i) + N(i);

y_1 = (HR_1(i)^2 + HI_1(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);
y_2 = (HR_2(i)^2 + HI_2(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);
y_3 = (HR_3(i)^2 + HI_3(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);



//demodulation steps
  
y_demod_I_1 =  y_1*(cos(2*%pi*fc*n*Ts))';
y_demod_Q_1 =  y_1*(sin(2*%pi*fc*n*Ts))';

y_demod_I_2 =  y_2*(cos(2*%pi*fc*n*Ts))';
y_demod_Q_2 =  y_2*(sin(2*%pi*fc*n*Ts))';

y_demod_I_3 =  y_3*(cos(2*%pi*fc*n*Ts))';
y_demod_Q_3 =  y_3*(sin(2*%pi*fc*n*Ts))';

I_demod_1(i) = sign(y_demod_I_1)
Q_demod_1(i) = sign(y_demod_Q_1)

//X_1= complex(I_demod_1(i),Q_demod_1(i));



I_demod_2(i) = sign(y_demod_I_2)
Q_demod_2(i) = sign(y_demod_Q_2);

///X_2= complex(I_demod_2(i),Q_demod_2(i));

I_demod_3(i) = sign(y_demod_I_3)
Q_demod_3(i) = sign(y_demod_Q_3)


no_of_1_I = nnz(members([I_demod_1(i),I_demod_2(i),I_demod_3(i)],1));

if(no_of_1_I >= 2) I_demod(i) = 1;
else I_demod(i) = -1;   
end

no_of_1_Q = nnz(members([Q_demod_1(i),Q_demod_2(i),Q_demod_3(i)],1));

if(no_of_1_Q >= 2) Q_demod(i) = 1;
else Q_demod(i) = -1
  
end
if (I(i) ~= I_demod(i)) or (Q(i) ~= Q_demod(i))
    num_error = num_error + 1;    
end

i = i+1;

end


Pe(j) = num_error/n_points;


j = j+1;
end

plot(10*log10((A.^2)/0.001) - noise_powerdBm,(Pe),'colo','green');
xtitle( 'Pe vs Tx SNR', 'Tx SNR', 'Pe');
legend(["No diversity","Antenna diversity", "Time diversity"]);

////////////////////////////////////////////////////////////////////////

