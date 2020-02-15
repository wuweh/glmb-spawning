%=== signal 're-generation' using the same X

X = scenario.X;
L = scenario.L;
track_list = scenario.track_list;
tbirth = scenario.tBirth;
N_true = scenario.N_true;

Z= cell(K,1); %state and observation sets
Z_tracks = Z; % for troubleshooting
Lmiss = cell(K,1);


tmp_P_D_vec=cell(K,1);
for k=1:K
    if miss_on
        tmp_P_D_vec{k}= P_D*ones(length(track_list{k}),1);
        for idxtgt=1:length(track_list{k})
            if abs(tbirth(track_list{k}(idxtgt))-k) <= 0
                tmp_P_D_vec{k}(idxtgt)= 1;
            end
        end
    else tmp_P_D_vec{k}= ones(length(track_list{k}),1);
    end
        
end

for k=1:K
    if N_true(k)> 0
        idx= find( rand(N_true(k),1) <= tmp_P_D_vec{k} ); 
        if length(idx) ~= N_true(k)
            tidx = 1:length(track_list{k});
            lidx = setdiff(tidx,idx);
                Lmiss{k} = L{k}(:,lidx);
        end
        Z{k}= gen_observation(model,X{k}(:,idx));
    end
    Z_tracks{k} = Z{k};
    N_c= poissrnd(lambda_c);    %no. of clutter points
    C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(z_dim,N_c);  %clutter generation
    Z{k}= [ Z{k} C ];
end
    
