//system specifications
f_samp = 2000*10^6;
fc = 800*10^6;
Ts = 1/f_samp;
n = linspace(1,2000,2000);

n_points = 10000;
I = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
Q = 2*grand(n_points, 1, "bin", 1, 0.5) - 1;
A = 0.001

R = 100*10^(-3)
PLdB = 128.1 + 37.6*log10(R);
PL = 10^(-PLdB/20);


//fading
HR = grand(n_points, 1, "nor", 0, 1/sqrt(2));
HI = grand(n_points, 1, "nor", 0, 1/sqrt(2));

h = complex(HR,HI);


//coef of cos2pifct and sin2pifct
X = A*complex(I,Q);

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


//modulation step, samples of exp(2pifct) multipled with X= I +jQ and real part taken
y = (HR(i)^2 + HI(i)^2)*PL*(real(X(i))*cos(2*%pi*fc*n*Ts) + imag(X(i))*sin(2*%pi*fc*n*Ts)) + grand(1 ,length(n), "nor", 0, sigma);


//demodultion step: corelator reciever kind , Th = 0;

y_demod_I =  y*(cos(2*%pi*fc*n*Ts))';
y_demod_Q =  y*(sin(2*%pi*fc*n*Ts))';

I_demod(i) = sign(y_demod_I);
Q_demod(i) = sign(y_demod_Q);


//error count
if (I(i) ~= I_demod(i)) or (Q(i) ~= Q_demod(i))
    num_error = num_error + 1;    
end

i = i+1;
end


Pe = num_error/n_points;
