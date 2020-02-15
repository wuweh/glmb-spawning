function [assignments,costs]= mbestwrap_updt_gibbsamp_tempered(P0,P0_tempd,m)
    
n1 = size(P0,1);
n2 = size(P0,2);

assignments= zeros(m,n1);
costs= zeros(m,1)';

currsoln= n1+1:2*n1; %use all missed detections as initial solution
%[tempinit,~]= mofo(P0_tempd); currsoln= tempinit(:)'; %use sampling without replacement to generate initial solution
%[tempinit,~]= munkres(P0_tempd-min(min(P0_tempd))); currsoln= tempinit(:)'; %use slow matlab based munkres to generate initial solution
%[tempinit,~]= assignmentoptimal(P0_tempd-min(min(P0_tempd))); currsoln= tempinit(:)'; %use C based optimal assignment to generate initial solution

assignments(1,:)= currsoln;
costs(1)=sum(P0(sub2ind(size(P0),1:n1,currsoln)));

% %do in log domain
% for sol= 2:m
%     for var= 1:n1
%         tempsamp= P0_tempd(var,:); %grab row of costs for current association variable
%         tempsamp(currsoln([1:var-1,var+1:end]))= -log(0); %lock out current and previous iteration step assignments except for the one in question
%         tempsamp= -tempsamp; %restore positive weights
%         tempsamp(isfinite(tempsamp))= tempsamp(isfinite(tempsamp))-logsumexp(tempsamp(isfinite(tempsamp))); %normalize weights
%         currsoln(var)= resample(exp(tempsamp),1); %generate next sample from conditional
%     end
%     assignments(sol,:)= currsoln;
%     costs(sol)= sum(P0(sub2ind(size(P0),1:n1,currsoln)));
% end

%do in real domain
for sol= 2:m
    for var= 1:n1
        tempsamp= exp(-P0_tempd(var,:)); %grab row of costs for current association variable
        tempsamp(currsoln([1:var-1,var+1:end]))= 0; %lock out current and previous iteration step assignments except for the one in question
        idxold= find(tempsamp>0); tempsamp= tempsamp(idxold);
        [~,currsoln(var)]= histc(rand(1,1),[0;cumsum(tempsamp(:))/sum(tempsamp)]); 
        %currsoln(var)= resample(tempsamp/sum(tempsamp),1);
        currsoln(var)= idxold(currsoln(var));
    end
    assignments(sol,:)= currsoln;
    costs(sol)= sum(P0(sub2ind(size(P0),1:n1,currsoln)));
end

[C,I,~]= unique(assignments,'rows');assignments= C;
costs= costs(I);
