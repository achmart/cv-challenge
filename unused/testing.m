a = [1.2; 3.3; 4.4];
b = [2.2; 1.3; 0.4];


tic
for i=1:100000
    c = kron(a,b);
end
toc

tic
for i = 1:100000
    c = [a(1)*b; a(2)*b; a(3)*b];
end
toc