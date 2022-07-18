%%
function X = Bi_Measure(A,B)
    Sc = [];
    for i = 1:size(A,1)
        ma = sum(A(i,:),2);%
        mb = sum(B(i,:),2);
        Sc(i) = ma/mb; 
    end  
    X = sum(Sc,2)/size(A,1);

end


        