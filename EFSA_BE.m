function [old_parcorrf_BE,delet_edge_BE,record_edges_BE,spareBE]=EFSA_BE(old_var,N,q,cutoff)

global graph_form

%% This code is modified based on the Junbai Wang's Mgraph algorithm. 
%% This code is used to predict the regulations from non-linear terms to genes.


initial_model=eye(q,q);
k_sets=1:q;
full_model=triu(ones(q,q));
[set_a1,set_a2]=find(initial_model);
set_a=[set_a1,set_a2];
ln_set_a=size(set_a,1);



%rest of edges


for i=1:11
    for j=1:11
        if full_model(i,j)==1
            full_model(i,j)=0;
        end
    end
end




for i=12:77
    for j=12:77
        if full_model(i,j)==1
            full_model(i,j)=0;
        end
    end
end



[set_b1,set_b2]=find(full_model);
set_b=[set_b1,set_b2];
ln_set_b=size(set_b,1);




need_test_sets=setdiff(set_b,set_a,'rows');


%compement set of set a


for i=1:ln_set_a
    comp_set_a{i}=setdiff(k_sets,set_a(i,:)); %complement of set a
end




idx=1; %index for recorded edges
temp_G=zeros(q,q);
new_G=temp_G;
test1=0;
aa=1;


while aa < 58
   
    
    fit_var=EFSA_GaussVariance(old_var,set_a,comp_set_a);
    
   

    %initial model device
    det_old=det(old_var);
    det_fit=det(fit_var);
    total_devs=abs(N*(log(det_old)-log(det_fit)));
    df=size(need_test_sets,1);
    total_p=Chisq(total_devs,df);
    

    
    
    new_model_set=set_a;
    index=1:df;
    devs=[];
    chtest_p=[];
    

     
     for i=1:df
        new_model_set(ln_set_a+1,:)=need_test_sets(i,:);
        
        
        
        new_comp_set = cell(1,66);
        
        for j=1:size(new_model_set,1)
            new_comp_set{j}=setdiff(k_sets,new_model_set(j,:));
        end
        
   
        
        new_fit_var=EFSA_GaussVariance(old_var,new_model_set,new_comp_set);
        

        
        det_old=det(fit_var);
        det_fit=det(new_fit_var);
        devs(i)=abs(N*(log(det_old)-log(det_fit)));
        chtest_p(i)=Chisq(devs(i),1);
        
      
        
    end 
    %devs
 
    chtest_p(find(chtest_p<cutoff));

   
    
    %add test decomposable 
    sortedDev=[];
    sorted_need_test_sets=[];
   sorted_chtest_p=[];
    
    [sortedDev, idx_sortDev]=sort(devs);
    sorted_need_test_sets=need_test_sets(idx_sortDev,:);
    sorted_chtest_p=chtest_p(idx_sortDev);
    len_of_dev=length(sortedDev);
    
   
    for i=1:len_of_dev
         temp_G=new_G;
         temp_add_edge=sorted_need_test_sets(end-i+1,:);
         temp_G(temp_add_edge(1),temp_add_edge(2))=1;
         temp_G(temp_add_edge(2),temp_add_edge(1))=1;
         if graph_form==1
             cond1=(isDecomposableG(temp_G) & sorted_chtest_p(end-i+1)<cutoff);
         else
             cond1=sorted_chtest_p(end-i+1)<cutoff;
         end
         if cond1
             break;
         end
    end
     
  
    
    
    if cond1
        %add it to the model
        %add_edge=need_test_sets(idx_maxDev,:);
        add_edge=temp_add_edge;
        new_G(temp_add_edge(1),temp_add_edge(2))=1;
        new_G(temp_add_edge(2),temp_add_edge(1))=1;
        temp_G=new_G;
        %end add in
        
        set_a(end+1,:)=add_edge;
        comp_set_a{end+1}=setdiff(k_sets,set_a(end,:));
        old_need_test_sets=need_test_sets;
        need_test_sets=setdiff(old_need_test_sets,add_edge,'rows');
        ln_set_a=size(set_a,1);
        record_old_var=fit_var;
        delet_edge_BE= size(need_test_sets,1)+1;
        record_edges_BE{idx}=num2str(add_edge);
        record_edges_BE'
        idx=idx+1;
        %temp_newModel=new_model  ;
    else
        record_old_var=fit_var;
        delet_edge_BE= size(need_test_sets,1);

        break;
    end
    
aa=aa+1;    
  
end 


old_corrf_BE=record_old_var./sqrt(diag(record_old_var)*diag(record_old_var)');
inv_old_var_BE=pinv(record_old_var);
old_parcorrf_BE=-inv_old_var_BE./sqrt(abs(diag(inv_old_var_BE)*diag(inv_old_var_BE)'));
spareBE=((abs(old_parcorrf_BE)>0.000000001 & abs(old_parcorrf_BE)~=1));
spareBE = double(spareBE);
last_edge_BE = str2num(record_edges_BE{57});
Last_a = last_edge_BE(1);
Last_b = last_edge_BE(2);
if spareBE(Last_a,Last_b) == 0
    spareBE(Last_a,Last_b) = 1;
end

if spareBE(Last_b,Last_a) == 0
    spareBE(Last_b,Last_a) = 1;
end


