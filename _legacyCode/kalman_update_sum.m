function [w_update,m_update,P_update] = kalman_update_sum(z,wz,H,b,R,wp,m,P)

%performs the Kalman update step for sum of gaussian likelihood and preds
%inputs:    wz    - vector of weights of gaussian components in transition
%           z     - vector array of observation points z
%           H     - matrix array of linear transformation matrices H in likelihood
%           b     - vector array of offset vectors b in likelihood
%           R     - matrix array of covariance matrices Q in likelihood
%           wp    - vector of weights of gaussian components in predicted density
%           m     - vector array of means of components in predicted density
%           P     - matrix array of covariance matrices P in predicted density
%         
%outputs:   w_update - vector of weights of components in updated density
%           m_update - vector array of means of components in updated density
%           P_update - matrix array of covariance matrices in updated density

% "vectorization: predicted mixture times individual likelihood component, then repeated by individual likelihood component"

x_dim = size(m,1); z_dim = size(z,1);

zlength = length(wz); 
plength = length(wp);
numterms = zlength*plength;

wp= wp(:);
wz= wz(:);


%FOR LOOP - no caching (std)
w_update = zeros(numterms,1);
m_update = zeros(x_dim,numterms);
P_update = zeros(x_dim,x_dim,numterms);

for idxz=1:zlength
    for idxp=1:plength

        [logqz_temp,m_temp,P_temp] = kalman_update_s(z,H(:,:,idxz),b(:,idxz),R(:,:,idxz),m(:,idxp),P(:,:,idxp));
        idx = (idxz-1)*plength+idxp;
        w_update(idx) = wz(idxz)*wp(idxp)*exp(logqz_temp);
        m_update(:,idx) = m_temp;
        P_update(:,:,idx) = P_temp;

    end
end

% %FOR LOOP - with 1 step previous caching (fast for repeated measurement updates on same prediction and likelihood)
% w_update = zeros(numterms,1);
% m_update = zeros(x_dim,numterms);
% P_update = zeros(x_dim,x_dim,numterms);
% 
% persistent mu_temp S_temp det_S_temp iS_temp K_temp; persistent old_wz old_H old_b old_R old_wp old_m old_P;
% 
% if ~isequal(old_wz,wz) || ~isequal(old_H,H) || ~isequal(old_b,b) || ~isequal(old_R,R) || ~isequal(old_wp,wp) || ~isequal(old_m,m) || ~isequal(old_P,P)
%     mu_temp= zeros(z_dim,zlength,plength); S_temp= zeros(z_dim,z_dim,zlength,plength); det_S_temp=zeros(zlength,plength); iS_temp= zeros(z_dim,z_dim,zlength,plength); K_temp= zeros(x_dim,z_dim,zlength,plength);
%     for idxz=1:zlength
%         for idxp=1:plength
%             
%             mu_temp(:,idxz,idxp) = H(:,:,idxz)*m(:,idxp);
%             S_temp(:,:,idxz,idxp)  = R(:,:,idxz)+H(:,:,idxz)*P(:,:,idxp)*H(:,:,idxz)'; Vs= chol(S_temp(:,:,idxz,idxp)); det_S_temp(idxz,idxp)= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS_temp(:,:,idxz,idxp)= inv_sqrt_S*inv_sqrt_S';
%             K_temp(:,:,idxz,idxp)  = P(:,:,idxp)*H(:,:,idxz)'*iS_temp(:,:,idxz,idxp);
%             
%         end
%     end
%     old_wz= wz; old_H= H; old_b= b; old_R= R; old_wp= wp; old_m= m; old_P= P;
% end
% 
% for idxz=1:zlength
%     for idxp=1:plength
%                 
%         zshift = z - b(:,:,idxz); logqz_temp = -0.5*z_dim*log(2*pi) - 0.5*log(det_S_temp(idxz,idxp)) - 0.5*(zshift-mu_temp(:,idxz,idxp))'*iS_temp(:,:,idxz,idxp)*(zshift-mu_temp(:,idxz,idxp));
%         
%         idx = (idxz-1)*plength+idxp;
%         w_update(idx) = wz(idxz)*wp(idxp)*exp(logqz_temp);
%         m_update(:,idx) = m(:,idxp) + K_temp(:,:,idxz,idxp)*(zshift-mu_temp(:,idxz,idxp));
%         P_update(:,:,idx) = (eye(size(P(:,:,idxp)))-K_temp(:,:,idxz,idxp)*H(:,:,idxz))*P(:,:,idxp);
%         
%     end
% end
% 
% 
% %VECTORIZATION - memory intensive due to large arrays
% zidxs= repmat(1:zlength,[plength 1]);  zidxs= zidxs(:);
% 
% wp_temp= repmat(wp,[zlength 1]);
% wz_temp= wz(zidxs);
% 
% mp_temp= repmat(reshape(m,[x_dim 1 plength]),[1 1 zlength]);
% Hz_temp= H(:,:,zidxs);
% 
% Pp_temp= repmat(P,[1 1 zlength]);
% Rz_temp= R(:,:,zidxs);
% 
% 
% mu_temp= mtimesx(Hz_temp,mp_temp); mu_temp= reshape(mu_temp,[z_dim,numterms]);
% S_temp= Rz_temp+ mtimesx(mtimesx(Hz_temp,Pp_temp),Hz_temp,'T');
% det_S_temp= real(prod(ndfun('eig',S_temp))); 
% iS_temp= multinv(S_temp); 
% iS_temp= (iS_temp+mtimesx(iS_temp,'T',1))/2;
% K_temp= mtimesx(mtimesx(Pp_temp,Hz_temp,'T'),iS_temp);
% 
% zs_temp= repmat(z,[1 numterms])-b(:,zidxs);
% eta_temp= zs_temp-mu_temp; eta_temp=reshape(eta_temp,[z_dim 1 numterms]);
% logqz_temp= -0.5*z_dim*log(2*pi) - 0.5*log(det_S_temp(:)) - 0.5*reshape(mtimesx(mtimesx(eta_temp,'T',iS_temp),eta_temp),[numterms 1]);
% 
% w_update= wz_temp.*wp_temp.*exp(logqz_temp(:));
% m_update= mp_temp + mtimesx(K_temp,eta_temp); m_update= reshape(m_update,[x_dim,numterms]);
% P_update= mtimesx(repmat(eye(x_dim),[1 1 numterms])-mtimesx(K_temp,Hz_temp),Pp_temp);

