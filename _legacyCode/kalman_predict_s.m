function [m_predict,P_predict] = kalman_predict_s(F,Q,m,P)

%perform the Kalman prediction
%int[N(x;F*xi,Q)*N(xi;m,P);d(xi)] = N(x;m_predict,P_predict)
%                                   = N(x;F*m,Q+FPF') 
%input:  F,Q,m,P 
%output: m_predict, P_predict

m_predict = F*m;
P_predict = Q + F*P*F';