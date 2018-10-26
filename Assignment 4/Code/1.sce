clear;
power_received_dBm = -60;
power_received_dB = power_received_dBm - 30;

BPSK_amp = sqrt(2*10^(power_received_dB/10));

noise_power_dBm = -60;
noise_power_dB = noise_power_dBm - 30;
No = 10^(noise_power_dB/10)
sigma = sqrt(No/2);


//no of Mobile stations:
M = linspace(10,100,10);

i = 1;

message_size = 1000; 
Pe = zeros(10,1);


while i <=length(M)
  
  //foreach M station, generate a random BPSK message sequence and a wideband random sequence
    message_seq = BPSK_amp*(2*grand(M(i), message_size, "bin", 1, 0.5)-1);
    PN_seq = 2*grand(M(i), 512, "bin", 1, 0.5)-1;
    decoded_code_word = zeros(M(i),message_size);
    j = 1;
    num_error = 0;
    while j <= message_size //message_index
        k = 1;
        code_word_sum = zeros(512,1);
        while k <= M(i)// user index
            
            code_word_sum = code_word_sum + (message_seq(k,j))*PN_seq((k-1)*512+1:(k*512))
           
            k = k+1
         end   
            received_code_word   = code_word_sum + grand(512 ,1, "nor", 0, sigma);
     
     //decoding
     k = 1;
     
     while k <= M(i)// user index
            
            decoded_code_word(k,j) = (received_code_word)'*(PN_seq((k-1)*512+1:(k*512)))
           
           
         if(sign(decoded_code_word(k,j)) ~= sign(message_seq(k,j))) num_error = num_error + 1;
         end
         
         
         
            k = k+1
     end
     

     j = j+1
    end
    
    Pe(i) = num_error/(message_size*(M(i)))
  
    i = i+1;
end

plot2d(M,Pe,logflag = "nl");

xtitle( 'BER vs No. of MSs (M)', 'No. of MSs (M)', 'BER');

