function [outputArg1] = create_cell(coor,n_cell,r_c)
[N,~]=size(coor);
 
cell_list=cell(n_cell,n_cell); 
for i=1:N
    xx=ceil(coor(i,1)/r_c);
    yy=ceil(coor(i,2)/r_c);
    for k=-1:1
        for l=-1:1
            cell_list{mod(xx-1+k,n_cell)+1,mod(yy-1+l,n_cell)+1}(end+1)=i;
        end
    end
end
outputArg1 = cell_list;
end