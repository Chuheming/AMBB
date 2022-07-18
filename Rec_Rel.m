%% evaluation 
%%Relevance and Recovery
% A represnts output cluster; B represents implanted cluster
% Rel = Bi_measure(A,B); 
% Rec = Bi_measure(B,A);

function X = Rec_Rel(A,B)
    Sc = [];
    [br,bc] = size(B);
    [ar,ac] = size(A);
    unsc = bc;
    for i = 1:ar
        sc = 0;
        A1 = sort(A(i,:));
        for j = 1:br
            
            B1 = sort(B(j,:));
            sc1 = size(intersect(A1,B1),2);
            if sc<sc1
                sc = sc1;
                unsc = size(unique(union(A1,B1)),2);
            end
        end
        Sc(i) = sc/unsc;
    end
    X = sum(Sc,2)/ar;
   
end