//Defining system parameters
lambda = 1;
mu = 1/120;
block_prob = [0,0,0,0,0,0,0];
N = [80,90,100,110,120,130,140];


for num_of_N = 1:7
    block_count = 0;

//As the arrival is a Poisson process of parameter(lamda) the 
//interarrival times are exponential with parameter lamda
//Thus genrating the inter arrival times and call hold times
        
        inter_arr_time = grand(1, 20000, "exp",lambda);
        
//arrival time calculation
        arr_time = zeros(1,20000);

        for i = 1:20000
            arr_time(i) = sum(inter_arr_time(1:i));
        end

//service time calculation
        service_time = grand(1,20000, "exp",1/mu);

        empty_mat = zeros(1,N(num_of_N));
//channel busy array stores the times till which the channels are going to be busy
        channel_busy_time = zeros(1,N(num_of_N));

//we compare the arrival time of the current user for the all elements of the channel busy array
//and check whether the min element is less than or greater than zero. Greater than zero imples call blocking  
        for k = 1:20000
                temp_mat = channel_busy_time - arr_time(k);
                
                if(min(temp_mat) > 0)
                    block_count = block_count + 1;
                else
                    [value,index] = min(channel_busy_time);
                    channel_busy_time(index) = arr_time(k) + service_time(k);
                end
        end
    
block_prob(num_of_N) = block_count/(20000);   
end
theo_prob = [0,0,0,0,0,0,0];
rho = lambda/mu;
for i = 1:7
    temp_sum = 0;
    for t = 1:N(i)
        temp_sum = temp_sum + rho^t/factorial(t);
    end
    theo_prob(i) = (rho^N(i))/(factorial(N(i))*temp_sum);
end

scf(2);
clf(2); 
plot(N,block_prob, 'ro-'); 
plot(N,theo_prob, 'bs:');
xtitle( 'Blocking Probability vs No. of frequency channels', 'Blocking probability', 'CDF')
legend(["Simulated Pr", "Theoretical Pr"]);




