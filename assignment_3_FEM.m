clear all;
close all;

%values

E = 1; %E-module
L = 1; %rod-length
q = 1; %distributed load
Aends = [2 1]; %area in the start and end of the rod
n = 6; %number of elements
Le = L/n; %element length
m = n + 1; % number of nodes

%the area of the rod in each node
A = zeros(1,m);%preallocating memory for the A-vector
for i = 1:m
    A(i) = Aends(1) + (Aends(2)-Aends(1))*(i-1)*Le/L;
end

%the stiffness matrices for each element made in multidimensional arrays
%each matrix is one layer in the array, and the i-th layer can be found
%by calling k(:,:,i)
k = zeros(2,2,n); %preallocating memory
for i = 1:n
    k(1:2,1:2,i) = [1 -1; -1 1]*E*(A(i)+A(i+1))/(2*Le);
end

%finding the nodal load vectors
S0 = zeros(2,1,n);%preallocating memory
for i = 1:n
    S0(1:2,1:1,i) = -q*Le/2*[1;1];
end

%finding the connectivity matrices
a = zeros(2,m,n);%preallocating memory for m DOFs
for i = 1:n
    a(1:2,i:i+1,i) = [1 0; 0 1];
end

%finding the global stiffness matrix
K = zeros(m,m);%allocating memory for a stiffness matrix with m DOFs

for i = 1:n
    K(1:m,1:m) = a(:,:,i).'*k(:,:,i)*a(:,:,i) + K;
end

%finding the consistent load vector
R0 = zeros(m,1);
for i = 1:n
    R0(1:m,1) = a(:,:,i).'* S0(:,:,i) + R0;
end

%load vector
R(m,1) = q*L/2; %load on the right end
R = R-R0;

%including BCs --> the DOF at the left end is fixed
R(1) = [];
K(1,:) =[];
K(:,1)=[];

%solved R=Kr for the displacements r
%mldivide : x = A\B produces the particular solution of the linear equation
%A*x = B
r = K\R;

%finds displacement and axial forci in each element of the rod
x(1,:) = 0:0.001:Le;
xx = length(x);
u(1,:) = x/Le*r(1);
N(1,:) = E*Aends(2)/Le*r(1)*(2-x/L);
for i = 2:n
    x(i,:) = (i-1)*Le*ones(1,xx) + (0:0.001:Le);
    
    u(i,:) = ((1-x(1,:)/Le)*r(i-1)+x(1,:)/Le*r(i));
    N(i,:) = (E*Aends(2)/Le)*([-1 1]*[r(i-1) r(i)].')*(2-x(i,:)/L);
end

%finds analytical solution
xa = 0:0.001:L;
Na = q*(3*L/2-xa);
ua = q*L/(E*Aends(2))*(xa+L/2*log(1-xa/(2*L)));

%creates plots of both displacements in the same figure
figure(1)
plot(xa,ua) 
hold on
for i = 1:n
    plot(x(i,:),u(i,:),'r--')
end

xlabel('Rod axis [x/L]');
ylabel('Displacement[q0*L²/EA]');
legend('Analytical solution','FE solution');

figure(2)
plot(xa,Na)
hold on
for i = 1:n
    plot(x(i,:),N(i,:),'r--')
end

xlabel('Rod axis [x/L]');
ylabel('Force [q0*L]');
legend('Analytical solution', 'FE solution');
