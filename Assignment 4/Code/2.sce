clear;


function [inside] = point_in_polygon(xpol, ypol, xpoint, ypoint)
//*****************************************************************************
// function: inside_poly - this function returns a Yes (1), No(0)
// depending on if a point is inside a  polygon or not
//
// Inputs:
//    xpol, ypol are the coordinates of a simple closed polygon 
//    xpoint, ypoint are the coordinates of the point
//
// Outputs: 1 (inside polygon), 0 (not inside polygon)
//
//*****************************************************************************

// make them row vectors instead of column... consistency
if(size(xpol,2)==1) 
    xpol=xpol';
end
if(size(ypol,2)==1) 
    ypol=ypol';
end

//********************************************
// show the polygon
do_the_plot = 0

if(do_the_plot)
    scf(9999);
    // it is assumed that this is a closed polygon
    plot2d(xtemp,ytemp)
end
//********************************************

npol = size(xpol,'*')


inside = 0
j = npol; // j is the previous vertice 
i = 1
while  i <= npol

    if ((((ypol(i) <= ypoint) & (ypoint < ypol(j)))|((ypol(j) <= ypoint) & (ypoint < ypol(i)))) & ..
       (xpoint < ((xpol(j) - xpol(i))/(ypol(j) - ypol(i))) * (ypoint - ypol(i)) + xpol(i)))

          inside = 1;
    end
    i = i + 1;
    j = i - 1;
end

endfunction
//Power transmitted
Pt = 30;
gamma1 = 3.5; 
Rc = 500;
R = (sqrt(3)/2)*Rc;
//coordinated of the center of the interfering cells
u = [1,1,0,-1,-1,0];
v= [0,-1,-1,0,1,1];
//coordinates of ref cell
x_ref = [1/2,1,1/2,-1/2,-1,-1/2]*(Rc);
y_ref = [(sqrt(3)/2),0,-(sqrt(3)/2),-(sqrt(3)/2),0,(sqrt(3)/2)]*(Rc);
//plot(x_ref,y_ref,'o');


//cordinates of center of interfering cells
x_intf = (sqrt(3)/2)*u*(2*R);
y_intf = (1/2)*u*(2*R) + v*(2*R);


// generate 5000 random points in the square which circumscribes the reference cell
//square
y_square1 = -Rc;
y_square2 = Rc;
x_square1 = -Rc;
x_square2 = Rc;
//random coordinate generation
x_r = grand(5000, 1, "unf", x_square1, x_square2);
y_r = grand(5000, 1, "unf", y_square1, y_square2); 
i=1;
j=1;
required_points_x = zeros(5000,1);
required_points_y = zeros(5000,1);

x_ref_neg = (-1)*x_ref;
y_ref_neg = y_ref;
x_r_neg = (-1)*x_r;
y_r_neg = y_r;

//only taking those points which fall inside the referece hexagon
while i <= 5000
   
   IN1 = point_in_polygon(x_ref, y_ref, x_r(i), y_r(i));
   IN2 = point_in_polygon(x_ref_neg, y_ref_neg, x_r_neg(i), y_r_neg(i));
   
   dist = sqrt((x_r(i))^2 + (y_r(i))^2);   
   if(IN1 & IN2 & (dist > 10)) 
       required_points_x(j) = x_r(i);
       required_points_y(j) = y_r(i);
       j = j+1;
   end    
   i = i+1;
end

//no of Mobile stations:
M = linspace(10,100,10);


power_received_dBm = -60;
power_received_dB = power_received_dBm - 30;

BPSK_amp = sqrt(2*10^(power_received_dB/10));

noise_power_dBm = -60;
noise_power_dB = noise_power_dBm - 30;
No = 10^(noise_power_dB/10)
sigma = sqrt(No/2);

message_size = 1000;
 user_index =1;
while user_index <= length(M)
  
  

required_points_refx = required_points_x(1:M(user_index));
required_points_refy = required_points_y(1:M(user_index));

//plot(x_ref1,y_ref1);
//plot(required_points_refx,required_points_refy,'o')
k = 1;

while k<= 6
    
    required_points_intx((k-1)*M(user_index) + 1: (k)*M(user_index)) = x_intf(k) + required_points_x(k*M(user_index) + 1: (k+1)*M(user_index));
    required_points_inty((k-1)*M(user_index) + 1: (k)*M(user_index)) = y_intf(k) + required_points_y(k*M(user_index) + 1: (k+1)*M(user_index));
    //plot(x_ref1 + x_intf(k) ,y_ref1 + y_intf(k));
    //plot((required_points_intx((k-1)*M(user_index) + 1: (k)*M(user_index))),(required_points_inty((k-1)*M(user_index) + 1: (k)*M(user_index))) ,'o');
    
    k = k+1;
end
  
  
  
  //foreach M station, generate a random BPSK message sequence and a wideband random sequence
    message_seq = BPSK_amp*(2*grand(M(user_index), message_size, "bin", 1, 0.5)-1);
    PN_seq = 2*grand(M(user_index), 512, "bin", 1, 0.5)-1;
    
    PN_seq_inter = 2*grand(6*M(user_index), 512, "bin", 1, 0.5)-1;
    message_seq_inter = (2*grand(6*M(user_index), message_size, "bin", 1, 0.5)-1);
    
    decoded_code_word = zeros(M(user_index),message_size);
    num_error = 0;
    j = 1;
    
    while j <= message_size //message_index
        k = 1;
        code_word_sum = zeros(512,1);
        
        //adding the signals from the own (reference cell)
        while k <= M(user_index)// user index
            
            code_word_sum = code_word_sum + (message_seq(k,j))*PN_seq((k-1)*512+1:(k*512))
           
            k = k+1
         end   
         
         //adding the signals due to the neighbouring cells
         
         sigma_x = 4;
         b = 1/sqrt(2);
         
         neighbouring_cell_index = 1;
         
         while neighbouring_cell_index <= 6*M(user_index)
             
            d_own = min(sqrt((required_points_intx(neighbouring_cell_index) - x_intf(1))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(1))^2),sqrt((required_points_intx(neighbouring_cell_index) - x_intf(2))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(2))^2),sqrt((required_points_intx(neighbouring_cell_index) - x_intf(3))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(3))^2),sqrt((required_points_intx(neighbouring_cell_index) - x_intf(4))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(4))^2),sqrt((required_points_intx(neighbouring_cell_index) - x_intf(5))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(5))^2),sqrt((required_points_intx(neighbouring_cell_index) - x_intf(6))^2 + (required_points_inty(neighbouring_cell_index) - y_intf(6))^2));
            
            d_ref = sqrt((required_points_intx(neighbouring_cell_index))^2 + (required_points_inty(neighbouring_cell_index))^2);
             phi_10 = b*grand(1,1,"nor",0,sigma_x);
             //P_inter = (10^(power_received_dB/10))*((d_own)^(3.5))*(10^(phi_10)/10)/((d_ref)^(3.5));
             P_inter = (10^(power_received_dB/10))*((d_own)^(3.5))/((d_ref)^(3.5));
             BPSK_amp_inter = sqrt(2*P_inter);
             
             //PN_seq_inter = 2*grand(512, 1, "bin", 1, 0.5)-1;
             //message_seq_inter = BPSK_amp_inter*(2*grand(1, 1, "bin", 1, 0.5)-1);
             
             code_word_sum = code_word_sum + (BPSK_amp_inter)*(message_seq_inter(neighbouring_cell_index,j))*PN_seq_inter((neighbouring_cell_index-1)*512+1:(neighbouring_cell_index*512));
             
             neighbouring_cell_index = neighbouring_cell_index +1;
         end
         
         
         
         
         
         
            received_code_word   = code_word_sum + grand(512 ,1, "nor", 0, sigma);
     
     //decoding
     k = 1;
     
     while k <= M(user_index)// user index
            
            decoded_code_word(k,j) = (received_code_word)'*(PN_seq((k-1)*512+1:(k*512)))
           
           
         if(sign(decoded_code_word(k,j)) ~= sign(message_seq(k,j))) num_error = num_error + 1;
         end
         
         
         
            k = k+1
     end
     

     j = j+1
    end
    
    Pe(user_index) = num_error/(message_size*(M(user_index)))
  
    user_index = user_index+1;
end

plot2d(M,Pe,logflag = "nl");
xtitle( 'BER vs No. of MSs (M)', 'No. of MSs (M)', 'BER');



