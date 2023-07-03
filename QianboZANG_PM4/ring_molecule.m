function [type_r, pos_r, vel_r, bond_r] = ring_molecule()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

type_r = [];
pos_r = [];
vel_r = [];
for i = 1:10
    type_r = [type_r;[1,1,1,1,1,1,1,1,1,1].'];
    
    tmp = 0.1*rand(10,2)+randi(10)+1;

    pos_r = [pos_r;tmp];
    vel_r = [vel_r;zeros(10,2)];
end

bond_r = []; % list of bonds
for i = 0:9
    bond_r = [bond_r; [10*i+1, 10*i+2]; [10*i+2, 10*i+3]; [10*i+3, 10*i+4]; [10*i+4, 10*i+5]; [10*i+5, 10*i+6]; [10*i+6, 10*i+7]; [10*i+7, 10*i+8]; [10*i+8, 10*i+9]; [10*i+9, 10*i+10]; [10*i+10, 10*i+1]];
end

end