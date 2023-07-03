clear;
clc;

% Constants
Nf = 900; 
L = 15; % box size
r_c = 1;
n_cell = ceil(L/r_c);
berendsen = 1;
T_b = 2;
k_s = 100; % spring constant
r_s = 0.3; % equilibrium bond length

dt = 0.01; % time step size
steps = 5000; % number of steps
aij = [50 25 25 200; 25 1 300 200; 25 300 25 200; 200 200 200 0]; 
gamma = 4.5; % coefficient for dissipative force
sigma = 1; % coefficient for random force

% Particle position (randomly initialized inside the box)
type_f = 3*ones(Nf,1);
pos_f = [(L-1)*rand(Nf,1)+0.5,L*rand(Nf,1)];
vel_f = 5*randn(Nf,2);

Nw = 1100;
type_w = 4*ones(200,1);
coor_w = [0.5*rand(100,1),L*rand(100,1);0.5*rand(100,1)+14.5,L*rand(100,1)];
vel_w = [zeros(100,1),5*ones(100,1);zeros(100,1),-5*ones(100,1)];

% chain mole
[type_l, pos_l, vel_l, bond_l] = ring_molecule();

% all particles
type = [type_f;type_w;type_l];
position = [pos_f;coor_w;pos_l];
velocity = [vel_f;vel_w;vel_l];
[N,~]=size(position);

bond_l = Nw + bond_l;

M = zeros(steps,2);
E_k=zeros(steps,1); 
T_c=zeros(steps,1); 

% link molecular

eps=randn(N,N);
% Make the matrix symmetric
epsilon=(eps+eps')/2;


for t = 1:steps
    Fc = zeros(N,2); % conservative force
    Fd = zeros(N,2); % dissipative force
    Fr = zeros(N,2); % random force
    Fs = zeros(N,2);

    cell_list=create_cell(position,n_cell,r_c);
    
    % Calculate forces
    for i = 1:N
        xx=ceil(position(i,1)/r_c);
        yy=ceil(position(i,2)/r_c);
        neighbor=cell_list{xx,yy};

        for j_idx=1:length(neighbor)
            j =  neighbor(j_idx);
            if j > i 

                a = aij(type(i),type(j));
        
                d_x = position(i,1)-position(j,1);
                d_y = position(i,2)-position(j,2);
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
                rij = [d_x,d_y];
                r = sqrt(d_x^2+d_y^2);
                
                % Conservative force
                if r < r_c
                    Fcij = a*(1-(r/r_c))*(rij/r);
                else
                    Fcij = [0,0];
                end
                Fc(i,:) = Fc(i,:) + Fcij;
                Fc(j,:) = Fc(j,:) - Fcij;
                
                % Dissipative and random forces
                if r < r_c
                    vij = velocity(i,:) - velocity(j,:);
                    wr = 1-r/r_c;
                    wd = wr^2;
                    Fdij = -gamma*wd*(dot(vij, rij)/r)*(rij/r);
%                     Frij = [(sigma*wr*epsilon(i,j)*d_x)/sqrt(dt), (sigma*wr*epsilon(i,j)*d_y)/sqrt(dt)];
                    Frij = sigma*sqrt(wr/dt)*randn(1,2); 
                else
                    Fdij = [0,0];
                    Frij = [0,0];
                end
                Fd(i,:) = Fd(i,:) + Fdij;
                Fd(j,:) = Fd(j,:) - Fdij;
                Fr(i,:) = Fr(i,:) + Frij;
                Fr(j,:) = Fr(j,:) - Frij;
                
            end
        end

        % bond force
        if i>Nw 
            for j = Nw:N
                if any(ismember(bond_l, [i, j], 'rows'))
                    d_x = position(i,1)-position(j,1);
                    d_y = position(i,2)-position(j,2);
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
         
                    r = sqrt(d_x^2+d_y^2);
                
                    Fsij = [k_s*(1-r/r_s)*d_x,k_s*(1-r/r_s)*d_y];
                    Fs(i,:) = Fs(i,:) + Fsij;
                    Fs(j,:) = Fs(j,:) - Fsij;

                end
            end          
        end
    end
    % reset the force of wall || add F-body
    for w = Nf+1: 1000
        Fr(w,:) = [0, 0.3];
        Fd(w,:) = [0, 0];
        Fc(w,:) = [0, 0];
    end

    for w = 1001: 1100
        Fr(w,:) = [0, -0.3];
        Fd(w,:) = [0, 0];
        Fc(w,:) = [0, 0];
    end


    % Update velocities
    velocity = velocity + dt*(Fc + Fd + Fr + Fs);

    % Update positions with periodic boundary conditions
    position = position + dt*velocity;
    position = position - L*floor(position/L);

    M(t,:)=sum(velocity);
    E_k(t)=sum(sum(velocity.*velocity))/2;
    T_c(t)=E_k(t)*(2/3)/N;

    if berendsen    % Berendsen thermostat
        velocity = velocity * sqrt(1+0.0025*(T_b/T_c(t)-1));
    end


    
    % Visualization
    if mod(t,10) == 0
        scatter(position(1:Nf,1), position(1:Nf,2), 'filled', 'red'); % plot particles
        hold on
        scatter(position(Nf+1:Nw,1), position(Nf+1:Nw,2),50, 'black', 'filled');
        hold on
        scatter(position(Nw+1:end,1), position(Nw+1:end,2), 'blue', 'filled');
        hold off
        axis([0 L 0 L]); % size of the box
        title(['Momentum: ', num2str(M(t,:))]);
        drawnow;
    end
end

% plot(1:steps, T_c);
% Me = M(:,1) + M(:,2);
% plot(1:steps, Me)
