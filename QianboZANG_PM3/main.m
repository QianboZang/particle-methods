L=30;           
N=100;           
T=0.1;             
r_c=2.5; 
berendsen=true;
T_b=1;
n_cell=ceil(L/r_c);
d_t=0.001;
     
num_epochs=5000;       
E_p=zeros(num_epochs,1); 
E_k=zeros(num_epochs,1);
E_total=zeros(num_epochs,1); 
T_c=zeros(num_epochs,1);   
M=zeros(num_epochs,2);    

p_coor=zeros(N,2);  
p_v=zeros(N,2);     
p_a=zeros(N,2);  
p_a_new=zeros(N,2);
% create a shuffler to make sure particle not overlaped.
shuffler=randperm(L^2);
for i=1:N
    t=shuffler(i)-1;
    p_coor(i,:)=[mod(t,L),fix(t/L)]; 
    p_v(i,:) = [rand()-0.5,rand()-0.5];
end
% make sure the started momentum is 0
ini_m=sum(p_v,1)/N;
p_v=p_v-ini_m;  

for e=1:num_epochs
    % update the coordanates of particles 
    p_coor=p_coor+p_v*d_t+p_a*d_t*d_t/2;  
    % periodic boundary conditions
    p_coor=mod(p_coor,L); 
    cell_list=create_cell(p_coor,n_cell,N,r_c);

    % compute potential and force/acceleration
    p_a=p_a_new;
    p_a_new=zeros(N,2);
    for i=1:N-1
        xx=ceil(p_coor(i,1)/r_c);
        yy=ceil(p_coor(i,2)/r_c);
        neighbor=cell_list{xx,yy};
        for j=1:length(neighbor)
            j_value=neighbor(j);
            if j_value>i    % make sure non-repeated computing
                [dist,d_x,d_y]=get_distance(p_coor,L,i,j_value);         
                if dist<=r_c
                    f=-48*dist^(-14)+24*dist^(-8);
                    p_a_new(i,:)= [p_a_new(i,1)-f*d_x/dist,p_a_new(i,2)-f*d_y/dist];
                    p_a_new(j_value,:)=[p_a_new(j_value,1)+f*d_x/dist,p_a_new(j_value,2)+f*d_y/dist];
                    E_p(e)=E_p(e)+4*dist^(-12)-4*dist^(-6)-4*r_c^(-12)+4*r_c^(-6);
                end
            end
        end
    end
    % update other values
    p_v=p_v+(p_a+p_a_new)*d_t/2; 
    M(e,:)=sum(p_v);
    E_k(e)=sum(sum(p_v.*p_v))/2;
    T_c(e)=E_k(e)*(2/3)/N;
    E_total(e)=E_p(e)+E_k(e);

    if berendsen    % Berendsen thermostat
        p_v=p_v*sqrt(1+0.0025*(T_b/T_c(e)-1));
    end

    if mod(e,10)==0
        pause(0.01) 
        scatter(p_coor(:,1),p_coor(:,2),20,'black')
    end 
end



