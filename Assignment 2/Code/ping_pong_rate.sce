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
stacksize('max');


//system parameters
Rc = 250;
R = (sqrt(3)/2)*Rc;

//coordinates of ref cell
x_ref = [(sqrt(3)/2),(sqrt(3)/2),0,-(sqrt(3)/2),-(sqrt(3)/2),0]*(Rc);
y_ref = [1/2,-1/2,-1,-1/2,1/2,1]*(Rc);
x_ref1 = [(sqrt(3)/2),(sqrt(3)/2),0,-(sqrt(3)/2),-(sqrt(3)/2),0,(sqrt(3)/2)]*(Rc);
y_ref1 = [1/2,-1/2,-1,-1/2,1/2,1,1/2]*(Rc)
//plot(x_ref1,y_ref1);
y_square1 = -Rc;
y_square2 = Rc;
x_square1 = -Rc;
x_square2 = Rc;


//random coordinate generation
x_r = grand(5000, 1, "unf", x_square1, x_square2);
y_r = grand(5000, 1, "unf", y_square1, y_square2); 
i=1;
j=1;
required_points_x = zeros(100,1);
required_points_y = zeros(100,1);

x_ref_neg = (-1)*x_ref;
y_ref_neg = y_ref;
x_r_neg = (-1)*x_r;
y_r_neg = y_r;

//The below generates 100 random points for the MS position in the reference cell. 
//only taking those points which fall inside the referece hexagon
while i <= 100
   
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
required_points_x = required_points_x(1:j-1);
required_points_y = required_points_y(1:j-1);

//plot(required_points_x,required_points_y,'o');

//bounding circle;
r_bounding = 3.12*Rc
a = linspace(0, 2*%pi, 100);
x_b = (r_bounding)*cos(a);
y_b = (r_bounding)*sin(a);
//plot(x_b, y_b);

//neighbouring BS center coordinates:
R = (sqrt(3)/2)*Rc;
ang_1 = [60,120,180,240,300,0];
ang_2 = [60,90,120,150,180,210,240,270,300,330,0];
ang_3 = [90,150,210,270,330,30];
x_neigh_1 = 2*R*cos(ang_1*%pi/180);
y_neigh_1 = 2*R*sin(ang_1*%pi/180);

x_neigh_2 = 4*R*cos(ang_1*%pi/180);
y_neigh_2 = 4*R*sin(ang_1*%pi/180);

x_neigh_3 = 3*Rc*cos(ang_3*%pi/180);
y_neigh_3 = 3*Rc*sin(ang_3*%pi/180);

x_neigh = [0,x_neigh_1,x_neigh_2,x_neigh_3];
y_neigh = [0,y_neigh_1,y_neigh_2,y_neigh_3];


//velocity (3km/hr = 5/6 m/s, 30km/hr = 50/6, 120km/hr= 200/6 )/////////////////////////////////////////////////////////
v = [5/6,50/6,200/6];
z = 3; //loop variable to iterate through velocity
H = [0,1,2,3,4];
// Direction theta (in degrees)
theta = grand(1, 1, "unf", 0, 360);

//array to keep count of the no. of ping pongs at different velocities
num_pingpong1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; 
num_pingpong2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
num_pingpong3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; 

//start of outer loop on velcity  
while z >=1
  
    //RSS measurement data;
    Pt = 30;
    h = 1; // loop variable to keep track of H
        
    // start of inner loop on H value
    while h <=5
        RSS_storage = hypermat([500000 19 1 1]);
        
        // initial starting point of the MS
        x_o = required_points_x(1);
        y_o = required_points_y(1);
        
        
        
        
        // movement coordinates
        num_steps = 100000;
        num_handouts = 0;
        
        x = zeros(num_steps,1); // initialization of all coordinates to (0,0)
        y = zeros(num_steps,1);
        
        x_handover = zeros(200,1);
        y_handover = zeros(200,1);
        
        x(1) = x_o;
        y(1) = y_o;
        t_step = 50*10^(-3);
        
        i = 2;
        j = 1;
        theta_initial = theta;
        current_BS = 1;
        m = 1; n = 1; p = 1;
        RSS_sum = zeros(1,19);
        
        prev_hand_x = 0;
        prev_hand_y = 0;
        prev_hand_t = 0;
        prev_BS = 0;
        n_ping_pong = 0;
        base_station_array = zeros(200,1);
        handover_times = zeros(200,1);
        random = 10*log10(grand(1,19, "exp",1))
        while num_handouts <=200//i <= num_steps
            
            x(i) = x(i-1) + v(z)*t_step*cos(theta*%pi/180);
            y(i) = y(i-1) + v(z)*t_step*sin(theta*%pi/180);
            
            //Reflection Logic
            //plot(x(i),y(i),'ro-');
            if((sqrt((x(i))^2 + (y(i))^2)) >= r_bounding)
              
                arrival_slope = y(i)/x(i);
                slope_initial = (y(i) - y_o)/(x(i) - x_o);
                angle = atan((slope_initial - arrival_slope)/(1+slope_initial*arrival_slope));
                if(j == 1)
                    d_o = sqrt((x(i) - x_o)^2 + (y(i) - y_o)^2);
                    d_s = sqrt((x_o)^2 + (y_o)^2);
                    alpha = acos(((d_o)^2 + (r_bounding)^2 - (d_s)^2)/(2*d_o*r_bounding))*180/%pi;
                    j = j+1;
                end
                a = theta;
                if(slope_initial > arrival_slope) 
                    theta = 180 - 2*alpha + theta;
                    a = 0;
                else
                    a = 1;   
                    theta = 180 + 2*alpha + theta;
                end
                x_o = x(i);
                y_o = y(i);        
            end
            
            //RSS measurements
            
            d = sqrt((x(i) - x_neigh).^2 + (y(i) - y_neigh).^2); // matrice of all the diatnces
            d_self = d(current_BS);
            
            RSS = Pt - (15.3 + 37.6*log10(d)) +10*log10(grand(1,19, "exp",1));
            RSS_storage(19*(i-1)+1 : 19*i) = RSS;
            
            
            if( i >10)
               
               RSS_sum = (RSS_storage(19*(i-1)+1 : 19*i) + RSS_storage(19*(i-2)+1 : 19*(i-1)) + RSS_storage(19*(i-3)+1 : 19*(i-2)) + RSS_storage(19*(i-4)+1 : 19*(i-3))+ RSS_storage(19*(i-5)+1 : 19*(i-4))+ RSS_storage(19*(i-6)+1 : 19*(i-5)) + RSS_storage(19*(i-7)+1 : 19*(i-6))+ RSS_storage(19*(i-8)+1 : 19*(i-7))+ RSS_storage(19*(i-9)+1 : 19*(i-8))+ RSS_storage(19*(i-10)+1 : 19*(i-9)))/10; 
            else  
            
            end    
            
            
            RSS_self = RSS_sum(current_BS);
            [val,index]=max(RSS_sum);
            
            // Handover Logic
            if (current_BS ~= index & (RSS_self + H(h) < val))
               
               current_BS = index;
               
                   
                  
               x_handover(m) = x(i);
               y_handover(m) = y(i);
               base_station_array(m) = index;
               handover_times(m) = t_step*i;
               //prev_hand_x = x(i);
               //prev_hand_y = y(i);
               //prev_hand_t = t_step*i;
              // plot(x(i),y(i),'o');
               m = m +1;
               num_handouts = num_handouts + 1;
              
               
            end    
                
            i = i+1;
            
        end
         k = 1;
         while k <=198
            
            if(base_station_array(k) == base_station_array(k+2))
            //if( base_station_array(k) == base_station_array(k+2))
                   if((handover_times(k+2) - handover_times(k+1) <= 1))
                     
                     if(z ==1)
                        num_pingpong1(h) = num_pingpong1(h) + 1;    
                     end    
                    if(z ==2)
                        num_pingpong2(h) = num_pingpong2(h) + 1;    
                     end
                   if(z ==3)
                        num_pingpong3(h) = num_pingpong3(h) + 1;    
                   end
                     
                   end
            end
            k = k+1;   
            
        end
         h = h+1;
     
     end
    
     
    
    
    
        
      z = z-1
end 

scf(2);
clf(2);
plot(H(1:5),num_pingpong1(1:5)/200,'ro-');
plot(H(1:5),num_pingpong2(1:5)/200,'bs:');
plot(H(1:5),num_pingpong3(1:5)/200);

//plot(v*18/5,num_handouts,'bs:');
xtitle( ' Ping Pong rates', 'Hysterisis', 'No. of handovers');
legend(["3km/hr","30km/hr", "120km/hr"]);


