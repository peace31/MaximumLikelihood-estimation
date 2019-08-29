function res = dummies(x)

L = length(x);
sorted_x = sort(x);

values_ind = [1; 1+find(sorted_x(1:L-1)~=sorted_x(2:L))];
values = sorted_x(values_ind);

D = length(values);
res = zeros(L,D);

for i = 1:D
    res(:,i) = (x==values(i));
end

end






