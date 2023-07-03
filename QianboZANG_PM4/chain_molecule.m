function [type_l, pos_l, vel_l, bond_l] = chain_molecule()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

type_l = [];
pos_l = [];
vel_l = [];
for i = 1:42
    type_l = [type_l;[1,1,2,2,2,2,2].'];
    
    tmp = 0.1*rand(7,2)+randi(10)+1;

    pos_l = [pos_l;tmp];
    vel_l = [vel_l;zeros(7,2)];
end

bond_l = []; % list of bonds
for i = 0:41
    bond_l = [bond_l; [7*i+1, 7*i+2]; [7*i+2, 7*i+3]; [7*i+3, 7*i+4]; [7*i+4, 7*i+5]; [7*i+5, 7*i+6]; [7*i+6, 7*i+7]];
end


end