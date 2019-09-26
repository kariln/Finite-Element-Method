%Assignment 5, rask 2 b)

clear all
close all
%number of elements, n
n = 6;

%creates the symbolic variables t0 and h
syms t0 h

%defines the total height H in terms of the n and h
H = h*n;

%creates R0 symbolically, with zero-values
R0 = sym(zeros(n+1,1));

%upgrading the values in R0 by iterating through the elements
for i = 1:n
    
    yc = H+h*(1-2*i);
    
    R0(i) = R0(i) - t0*h/(3*H^2)*(3*(H^2-yc^2)-h^2+2*h*yc);
    R0(i+1) = R0(i+1) -t0*h/(3*H^2)*(3*(H^2-yc^2)-h^2-2*h*yc);
end

Sum = sum(R0)
