%%
%%通过分析合成数据集，可知Match Score与可用数据占整体数据(1的总量)的百分比相等
clc
clear
time_t1 = cputime;
DDNum = 1;
Data = load ('Data_19_Engel.mat');
V = Data.in_X;
Genes = Data.Genes;
if size(Genes,2)>1
    Genes = Genes';
end
% V= load('text.txt');
X = getData(V);

% ori = load('Original.txt');

% X = Text_SVD(V);
bicluster = [];
score = 0;

num_bicluster = 0;%计算粗数

C = {}; %保存簇
cci = 1;

    sum_Xvalue = 0;
  
    m_x = 1;
    Match_score_Matric = [];
    Recovery = [];
    for row_threshould = 1500:1550
     
%     col_threshould = 3;


    bool = 1;
    N_T = X;
    cluster_num = 1;
    result_sum = 1;

    %获取行差异矩阵
    temp_N = getrow(N_T);
    %记录行标,即每行代表一个行簇
    RR = getnum_row(N_T,temp_N,row_threshould);
    
    r_k = [];
    max_col = [];
    i_row =1;
    m_y = 1;
    while (i_row<=size(RR,1))
        if size(RR,2)>1
        X_NT = [];
             if (RR(i_row,2)~= 0)
                %构造相似矩阵，将RR中的每一行的行簇元素存入X_NT中
                X_NT = Build_row(RR,N_T,i_row);
             end

        Max_col = get_max_col(X_NT);%获取列向量值的最大和

        col_threshould = size(X_NT,1)/2;
        
        [C_M,c_k] = getcol(X_NT,Max_col,col_threshould);%列向量聚类

        Boo = 1;
        if ~isempty(C_M)
        while Boo

            [C_M,r_k,c_k]=iter(C_M,RR(i_row,:),c_k);
            [p,q]=size(C_M);
            if sum(sum(C_M,1),2)==p*q    
                Boo = 0; 
                num_bicluster = num_bicluster + 1;
            end
        end
        end


        if ~isempty(C_M)&~isempty(r_k)&~isempty(c_k)
            r_k=r_k(1,:);
            X1_ori = [];

            for rki = 1:size(r_k,2)
                for k = 1:size(c_k,2)%获取聚类簇所在的行
                    X1_ori(rki,k) = V(r_k(rki),c_k(k));
                end
            end
           
            % 模拟数据集评估
%             result_Sc = Rec_Rel(X1_ori,ori);
%             if result_Sc>0
%                 Rec = Rec_Rel(ori,X1_ori);

%                 Match_score_Matric(m_x,m_y)=result_Sc;
%                 Recovery(m_x,m_y) = Rec;
%                 m_y = m_y + 1;
            % 真实数据集获取簇
                boo = test_cluste1(r_k,C);
                if boo == 1 && size(r_k,2)>=10
                    C{cci,1} = r_k;
                    C{cci,2} = c_k;
                    C{cci,3} = X1_ori;
                    cci = cci + 1;
                end
                
%             end
        end
        end
        i_row = i_row +1;

    end
%         Match_score_Matric = sort(Match_score_Matric,'Descend');
%         Recovery = sort(Recovery,'Descend');
    end

% Mean_Rele = mean(Match_score_Matric)
% Mean_Rec = mean(Recovery)

r = size(C,1);
for i = 1:r
    k = C{i,2};
    kk = size(k,2);
    t = 1;
    c = [];
    for j = 1:kk
        c{t,1} = Genes{k(1,j),1};
        t = t + 1;
    end
    C{i,4} = c;
end

filename = strcat('C.mat');
save(filename,'C');
time_t2 = cputime;
tmie_t = time_t2-time_t1

%%
%构造行差异矩阵
function  R_M = getrow(X_N)
    v_rm = [];
    for i= 1:size(X_N,1)
        for j = 1:size(X_N,1)
            v_rm(i,j) = sum(abs(X_N(i,:)-X_N(j,:)),2);
        end
    end
    R_M = v_rm;
end

%%
%记录行标
function Row_num = getnum_row(V,temp_N,row_threshould)
    
    i = 1;
    %Row = [];%保存行标
    
   Row = [];
   ri = 1;
    while(i<size(V,1))
        temp_row = 1;
        Row1 = [];
        Row1(1,temp_row) = i;%%记录种子下标
            for j = i+1:size(temp_N,2)    
               if temp_N(i,j) <= row_threshould
                   temp_row = temp_row +1;
                   Row1(1,temp_row) = j; % 记录所有符合阈值的行下标
               end
            end 
        i = i+1;
        boo = test_cluste(Row1,Row);
        if boo == 1
            nn = size(Row1,2);
            for nin = 1:nn
                Row(ri,nin) = Row1(1,nin);
                
            end
            ri = ri + 1;
        end
    end
    Row_num = Row;
end

%%
%获取聚类簇
function [C_M,c_k] = getcol(X_NT,max_score,threshold_col)%c_k返回的是列标号
%    while (max_score>1)
     c_k = [];
     i_t = 1;
     C_M1 = [];
     tmp_ck = [];
     boo_gc = 1;
     i = 1;
     while i <=  size(X_NT,2)&boo_gc     
        t = 1;
        if sum(X_NT(:,i),1) == max_score
            for j = 1 :size(X_NT,2)
            	if sum(abs(X_NT(:,i)-X_NT(:,j))) <= threshold_col 
                	c_k(i_t,t) = j;   %记录符合规则的列标
                	t = t+1;
                end
               	if size(c_k,2)>1
                  	boo_gc = 0;
                end         
            end
           	if size(c_k,2)>1&c_k(i_t,2)~=0
                tmp_ck = c_k(i_t,:);
               	C_M1 = getmatrix_col(tmp_ck,X_NT);
            end
           	if t>1
              	i_t = i_t +1;
            end
        end
     	i = i+1;
    end
    %得到最终双聚类簇
    if ~isempty(tmp_ck)
      	c_k = tmp_ck;
    else
      	c_k = [];
    end
    C_M = C_M1;
end
  
%%
 %获取聚类簇
function New_COL = getmatrix_col(c_k,X_NT)
    m = [];
    for i = 1:size(c_k,2)
        if c_k(i)~= 0
            m(:,i) = X_NT(:,c_k(i));
        end
    end
    New_COL = m;
end

%获取最大列向量和的值
function max_C = get_max_col(V)
max = 0;
    for i = 1:size(V,2)
        
        if max < sum(V(:,i),1)
            max = sum(V(:,i),1);
        end
       
    end
    max_C = max;
end

%%
function Res = Build_row(V,N_T,i_row)%输入V行差异矩阵，N_T矩阵，i_rowV中第几行
    X_NT = [];
    x_r = [];
    for q = 1:size(V(i_row,:),2)
        if V(i_row,q)~=0
            X_NT(q,:) = N_T(V(i_row,q),:);%X_NT为相似矩阵
            %x_r(q) = V(q);
        end
    end

    Res = X_NT;

end

%%
function [ite,Num_r,New_ck] = iter(C_M1,t_N1,c_k)%输入矩阵,t_N1表示行标矩阵,输出聚类簇，行标与列标
    R = [];
    boo_rn = 1;
    temp_N1 = getrow(C_M1); %行差异矩阵
    tr = 0;
    while boo_rn 
        Row_num1 = getnum_row(C_M1,temp_N1,tr);  %获取行标
        if size(Row_num1,2)<=1&&~isempty(Row_num1)%当只有一列时如何处理
            tr = tr +1;
        else
            Row_num1 = Row_num1(1,:);
            Row_num1(find(Row_num1==0)) = [];
            boo_rn = 0;
            f = size(Row_num1,2);
            FM= [];
            for fi = 1:f
               FM(fi,:) = C_M1(fi,:);
            end
            C_M1 = FM;
        end
    end
%构造列差异值矩阵获取列向量
    ms = get_max_col(C_M1) ; %最大值
    c_k2 = [];
    if ms>1
        [R,c_k1] = getcol(C_M1,ms,0);
        c_k_num = size(c_k1,2);
        
        if c_k_num>1
           for i = 1:c_k_num
                c_k2(i) = c_k(c_k1(i));
           end 
        end
%     else 
%         while i<=size(Row_num1,1)
%             R = Build_row(Row_num1,C_M1,i);
%             if sum(sum(cache_R,1),2)<sum(sum(R,1),2)
%                 cache_R = R;
%                 Julu_row = Row_num1(i,:);
%             end
%             i = i+1;
%         end
%         Num_r = getN_r(Julu_row,t_N1);
    end   
    t_N1(find(t_N1 == 0)) =[];
    Num_r = t_N1;

    New_ck = c_k2;
    ite = R;
end

%% 筛选
function New_c = test_cluste(r,R)
    rnum = size(r,2);
    clur = size(R,1);
    r = sort(r);
    boo = 0;
    difnum = 0;
    for i = 1:clur
        jb = 1;
        row = sort(R(i,:));
        k = size(row,2);
        for j = 1:rnum
            for jj = 1:k
                if row(1,jj) == r(1,j) %进度：检测不同的簇
                    break;  
                end
                if jj == k
                    jb = 0;
                end
            end
        end
        if jb == 0
            difnum = difnum + 1;
        end
    end
    if difnum == clur
        boo = 1;
    end
    New_c = boo;
end

%% test_cluste1(r_k,C)
function New_c = test_cluste1(r,R)
    rnum = size(r,2);
    clur = size(R,1);
    r = sort(r);
    boo = 0;
    difnum = 0;
    for i = 1:clur
        jb = 1;
        row = sort(R{i,1});
        k = size(row,2);
        for j = 1:rnum
            for jj = 1:k
                if row(1,jj) == r(1,j) %进度：检测不同的簇
                    break;  
                end
                if jj == k
                    jb = 0;
                end
            end
        end
        if jb == 0
            difnum = difnum + 1;
        end
    end
    if difnum == clur
        boo = 1;
    end
    New_c = boo;
end