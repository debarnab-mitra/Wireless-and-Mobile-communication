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
Rc = 1000;
R = (sqrt(3)/2)*Rc; 
//u,v coordinated of the center of the interfering cells
u = [2,3,1,-2,-3,-1];
v= [1,-2,-3,-1,2,3];
//coordinates of ref cell
x_ref = [1/2,1,1/2,-1/2,-1,-1/2]*(Rc);
y_ref = [(sqrt(3)/2),0,-(sqrt(3)/2),-(sqrt(3)/2),0,(sqrt(3)/2)]*(Rc);
//plot(x_ref,y_ref,'o');


//x, y cordinates of center of interfering cells
x_intf = (sqrt(3)/2)*u*(2*R);
y_intf = (1/2)*u*(2*R) + v*(2*R);

// generate 5000 random points in the square which circumscribes the reference cell
//square
y_square1 = -Rc;
y_square2 = Rc;
x_square1 = -Rc;
x_square2 = Rc;

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
required_points_x = required_points_x(1:j-1);
required_points_y = required_points_y(1:j-1);    
//plot(required_points_x,required_points_y,'o');
//take 1000 of those points which lie within the reference cell

i = 1;
SIR = zeros(j-1,1)
//calculate the SIR foreach of the considered points in the reference cell 
while i <= j-1

    
    dist = sqrt((required_points_x(i))^2 + (required_points_y(i))^2);
//intereference distances
    dist_1 = sqrt((required_points_x(i) - x_intf(1))^2 + (required_points_y(i) - y_intf(1))^2);
    dist_2 = sqrt((required_points_x(i) - x_intf(2))^2 + (required_points_y(i) - y_intf(2))^2);
    dist_3 = sqrt((required_points_x(i) - x_intf(3))^2 + (required_points_y(i) - y_intf(3))^2);
    dist_4 = sqrt((required_points_x(i) - x_intf(4))^2 + (required_points_y(i) - y_intf(4))^2);
    dist_5 = sqrt((required_points_x(i) - x_intf(5))^2 + (required_points_y(i) - y_intf(5))^2);
    dist_6 = sqrt((required_points_x(i) - x_intf(6))^2 + (required_points_y(i) - y_intf(6))^2);

//calculate SIR
    SIR(i) = ((1/(dist))^(gamma1))/((1/(dist_1))^(gamma1)+(1/(dist_2))^(gamma1)+(1/(dist_3))^(gamma1)+(1/(dist_4))^(gamma1)+(1/(dist_5))^(gamma1)+(1/(dist_6))^(gamma1));
    i = i+1;
end
// calculate the CDF of the SIR taken in dB scale
SIR_dB = 10*log10(SIR);
CDF = zeros(1000,1)
p = linspace(0,90,1000);
i = 2;
while i <=1000
    
    CDF(i) = (1/(j-1))*(nnz(dsearch(SIR_dB, [0,p(i)])))
    i = i+1
end

plot(p,CDF,'g');
xtitle( 'CDF for reuse 7', 'SIR', 'CDF')
