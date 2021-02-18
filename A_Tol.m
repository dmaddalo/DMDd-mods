function [eps1,eps11] = A_Tol(d,V)

[~,K]=size(V);

%% STEP 1: SVD of the original data

[~,Sigma,T]=svd(V,'econ');
sigmas=diag(Sigma);
Sigma = sparse(Sigma);
n=length(sigmas);

NormS=norm(sigmas,2); RMS1 = zeros(1,K);
for k=1:n
    RMS1(k) = norm(sigmas(k:n),2)/NormS;
end

figure
loglog(RMS1);
grid on; grid minor;
xlabel('mode'); ylabel('Reconstruction error');
title('SVD1 spectrum');

eps1 = input('First tolerance \n');
kk = find(RMS1 > eps1); kk = kk(end);

%% Spatial complexity: kk
('Spatial complexity')
kk

%% Create reduced snapshots matrix
hatT=Sigma(1:kk,1:kk)*T(:,1:kk)';
[N,~]=size(hatT);

%% Create the modified snapshot matrix
tildeT=zeros(d*N,K-d+1);
for ppp=1:d
    tildeT((ppp-1)*N+1:ppp*N,:)=hatT(:,ppp:ppp+K-d);
end

%% Dimension reduction
[~,Sigma1,~]=svd(tildeT,'econ');

clear tildeT V T

sigmas1 = diag(Sigma1);
n=length(sigmas1);

NormS=norm(sigmas1,2); RMS2 = zeros(1,n);
for k=1:n
    RMS2(k)=norm(sigmas1(k:n),2)/NormS;
end

figure
loglog(RMS2);
grid on; grid minor;
xlabel('mode'); ylabel('Reconstruction error');
title('SVD2 spectrum');

eps11 = input('Second tolerance \n');
kk1 = find(RMS1 > eps11); kk1 = kk1(end);

('Spatial dimension reduction')
kk1

 variables = whos; howbig = zeros(size(variables,1),1);
 for i = 1:size(variables,1)
    howbig(i) = variables(i).bytes;
 end
 howbigall = sum(howbig);
 
 fprintf(['So far you''ve used ',sprintf('%0.2f',howbigall/1e9),' GB of RAM \n']);
 if howbigall > 1e10
     fprintf('Beware!');
 end


end