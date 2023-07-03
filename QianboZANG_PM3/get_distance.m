function [outputArg1,outputArg2,outputArg3] = get_distance(p_coor,L,i,j_value)

d_x=p_coor(i,1)-p_coor(j_value,1);
d_y=p_coor(i,2)-p_coor(j_value,2);
if d_x>L/2 
    d_x=d_x-L;
elseif d_x<-L/2 
    d_x=L+d_x;
end
if d_y>L/2 
    d_y=d_y-L;
elseif d_y<-L/2 
    d_y=L+d_y;
end

outputArg1=sqrt(d_x^2+d_y^2);
outputArg2=d_x;
outputArg3=d_y;
end