function [logqz,m_update,P_update] = kalman_update_s(z,H,b,R,m,P)

%performs the Kalman update step
%N(z;Hx+b,R)*N(x;m,P) = qz*N(x,m_update,P_update)
%
%input:  z,H,b,R,m,P
%output: logqz,m_update,P_update

mu = H*m;
S  = R+H*P*H'; Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K  = P*H'*iS;

zshift = z - b;
z_dim = length(z);
logqz = -0.5*z_dim*log(2*pi) - 0.5*log(det_S) - 0.5*(zshift-mu)'*iS*(zshift-mu);
m_update = m + K*(zshift-mu);
P_update = (eye(size(P))-K*H)*P;

