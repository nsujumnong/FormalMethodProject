function [W, b] = generate_regression_mat_2(q1_data, q2_data, h)
syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;

l = size(q1_data, 1);
w = size(h,2);
W = zeros(l*2, w);
b = zeros(l*2, 1);
for i = 1:l
    W(2*(i-1)+1:2*(i-1)+2, :) = double(subs(h, {q1, q2, dq1, dq2, ddq1, ddq2},...
        {q1_data(i,1), q2_data(i,1),...
        q1_data(i,2), q2_data(i,2),...
        q1_data(i,3), q2_data(i,3)}));
    b(2*(i-1)+1:2*(i-1)+2, :) = [q1_data(i,4); q2_data(i,4)];
end
end

