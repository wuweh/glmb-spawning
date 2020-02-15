function [w_predict,m_predict,P_predict] = kalman_predict_sum(wt,F,Q,wp,m,P)

%performs the Kalman prediction for sum of gaussian transition and prior
%inputs:    wt    - vector of weights of gaussian components in transition
%           F     - matrix array of linear transformation matrices F in transition
%           Q     - matrix array of covariance matrices Q in transition
%           wp    - vector of weights of gaussian components in previous density
%           m     - vector array of means of components in previous density
%           P     - matrix array of covariance matrices P in previous density
%         
%outputs:   w_predict - vector of weights of components in predicted density
%           m_predict - vector array of means of components in predicted density
%           P_predict - matrix array of covariance matrices in predicted density
        
% "vectorization: posterior mixture times individual transition component, then repeated by individual transition component"

x_dim = size(m,1);

tlength = length(wt); 
plength = length(wp);
numterms = tlength*plength;

wp= wp(:);
wt= wt(:);

%FOR LOOP
w_predict = zeros(numterms,1);
m_predict = zeros(x_dim,numterms);
P_predict = zeros(x_dim,x_dim,numterms);

for idxt=1:tlength
    for idxp=1:plength
        
        [m_temp,P_temp] = kalman_predict_s(F(:,:,idxt),Q(:,:,idxt),m(:,idxp),P(:,:,idxp));
        idx = (idxt-1)*plength+idxp;
        w_predict(idx) = wt(idxt)*wp(idxp);
        m_predict(:,idx) = m_temp;
        P_predict(:,:,idx) = P_temp;
        
    end
end
% 
% %VECTORIZATION
% tidxs= repmat(1:tlength,[plength 1]);  tidxs= tidxs(:);
% 
% wp_temp= repmat(wp,[tlength 1]);
% wt_temp= wt(tidxs);
% 
% mp_temp= repmat(reshape(m,[x_dim 1 plength]),[1 1 tlength]);
% Ft_temp= F(:,:,tidxs);
% 
% Pp_temp= repmat(P,[1 1 tlength]);
% Qt_temp= Q(:,:,tidxs);
% 
% w_predict= wt_temp.*wp_temp;
% m_predict= mtimesx(Ft_temp,mp_temp); m_predict= reshape(m_predict,[x_dim,numterms]);
% P_predict= Qt_temp + mtimesx(mtimesx(Ft_temp,Pp_temp),Ft_temp,'T');


