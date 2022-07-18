%% binary matrix
function X = getData(V)
% Data = load ('D:\CaChe\Matlab_project\Data_Treutlin.mat');
    [m,n] = size(V);
    P = [];
    for i = 1:m
        xmax = max(V(i,:));
        xmin = min(V(i,:));
%         threshould = xmin+(xmax+xmin)/2;
        threshould = xmin+(xmax-xmin)/2;
        for j = 1:n
            if V(i,j) < threshould
                P(i,j) = 0;
            else
                P(i,j)=1;
            end
        end
    end
    X = P;
end
