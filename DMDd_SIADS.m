
 % This file is part of HODMD
 %
 % Copyright (c) 2017 S Le Clainche & J M Vega
 % All rights reserved.
 %
 % Redistribution and use in source and binary forms, with or without
 % modification, are permitted provided that the following conditions
 % are met:
 % 1. Redistributions of source code must retain the above copyright
 %    notice, this list of conditions and the following disclaimer.
 % 2. Redistributions in binary form must reproduce the above copyright
 %    notice, this list of conditions and the following disclaimer in the
 %    documentation and/or other materials provided with the distribution.
 %
 % THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 % ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 % ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 % OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 % HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 % LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 % OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 % SUCH DAMAGE.
 %
 % $FreeBSD$
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % FIND OUT HOW MUCH DISK SPACE IS BEING USED
 % 
%  variables = whos; howbig = zeros(size(variables,1),1);
%  for i = 1:size(variables,1)
%     howbig(i) = variables(i).bytes;
%  end
%  howbigall = sum(howbig);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [Vreconst,deltas,omegas,amplitude,modes] =DMDd_SIADS(d,V,Time,varepsilon1,varepsilon2,varepsilon)


%%%%%%%%%%%%%%%%%%%%%%%%%  DMD-d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function solves the HODMD algorithm presented in               %%%
%%% Le Clainche & Vega, SIAM J. on Appl. Dyn. Sys. 16(2):882-925, 2017  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %% INPUT: %%
%%% d: parameter of DMD-d (higher order Koopman assumption)
%%% V: snapshot matrix
%%% Time: vector time
%%% varepsilon1: first tolerance (SVD)
%%% varepsilon: second tolerance (DMD-d modes)
%%% %% OUTPUT: %%
%%% Vreconst: reconstruction of the snapshot matrix V
%%% deltas: growht rate of DMD modes
%%% omegas: frequency of DMD modes(angular frequency)
%%% amplitude: amplitude of DMD modes
%%% modes: DMD modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[J,K]=size(V);

%% STEP 1: SVD of the original data

[U,Sigma,W]=svd(V,'econ');
sigmas=diag(Sigma);
Sigma = sparse(Sigma);
n=length(sigmas);

NormS=norm(sigmas,2); RMS1 = zeros(1,K);
kk=0;
for k=1:n
    RMS1(k) = norm(sigmas(k:n),2)/NormS;
    if RMS1(k)>varepsilon1
        kk=kk+1;
    end
end

U=U(:,1:kk);

%% Spatial complexity: kk
('Spatial complexity')
kk

%% Create reduced snapshots matrix
hatV=Sigma(1:kk,1:kk)*W(:,1:kk)';
[N,~]=size(hatV);

%% Create the modified snapshot matrix
tilV=zeros(d*N,K-d+1);
for ppp=1:d
    tilV((ppp-1)*N+1:ppp*N,:)=hatV(:,ppp:ppp+K-d);
end

%% Dimension reduction
[tilU,tilSigma,tilW]=svd(tilV,'econ');

clear tilV

sigmas1=diag(tilSigma);
tilSigma = sparse(tilSigma);

Deltat=Time(2)-Time(1);
n=length(sigmas1);

NormS=norm(sigmas1,2); RRMSEE = zeros(1,n);
kk1=0;
for k=1:n
    RRMSEE(k)=norm(sigmas1(k:n),2)/NormS;
    if RRMSEE(k)>varepsilon2
        kk1=kk1+1;
    end
end

('Spatial dimension reduction')
kk1

tilU=tilU(:,1:kk1);
hattilV=tilSigma(1:kk1,1:kk1)*tilW(:,1:kk1)';

clear tilW

%% Reduced modified snapshot matrix
[hattilU,hattilSigma,hattilW]=svd(hattilV(:,1:end-1),'econ');

%% Reduced modified Koopman matrix
tilS=hattilV(:,2:end)*hattilW*inv(hattilSigma)*hattilU';
[tilQ,tilMM]=eig(tilS);
eigenvalues=diag(tilMM);

M=length(eigenvalues);
qq=log(eigenvalues);
deltas=real(qq)/Deltat;
omegas=imag(qq)/Deltat;

Q=tilU*tilQ;
Q=Q((d-1)*N+1:d*N,:);
[NN,MMM]=size(Q);

for m=1:MMM
    NormQ=Q(:,m);
    Q(:,m)= Q(:,m)/norm(NormQ(:),2);
end

%% Calculate amplitudes
% Mm=zeros(NN*K,M);
% Bb=zeros(NN*K,1);
% aa=eye(MMM);
% 
% clear U1 tilQ tilW hattilU hattilW tildeR hattilSigma hattilV
% 
% for k=1:K
%     Mm(1+(k-1)*NN:k*NN,:)=Q*aa;
%     aa=aa*tilMM;
%     Bb(1+(k-1)*NN:k*NN,1)=hatV(:,k);
% end
% 
% clear tildeMM hatT
% 
% CMm = Mm'*Mm;
% Trick = Mm'*Bb;
% 
% clear Mm Bb
% 
% [Vr,Sigmar2,~] = svd(CMm);
% 
% a = Vr*inv(Sigmar2)*Vr'*Trick;
% 
% % [Ur,Sigmar,Vr]=svd(Mm,'econ');
% % 
% % clear Mm
% % 
% % a=Vr*(Sigmar\(Ur'*Bb));
% 
% clear Vr Sigmar Ur Bb

%% Calculate amplitudes (lighter)

% Phi = tilU*hattilV; Phi = U*Phi((d-1)*N+1:d*N,2:end);
% [~,SigmaPhi,WPhi] = svd(Phi,'econ');

r = rank(tilS);
% KK = size(hattilW,1);

Vand = zeros(r,K);
zdmd = eigenvalues;

for i = 1:K
	Vand(:,i) = zdmd.^(i-1);
end

L = Q;
R = Vand;
G = hatV;
% G = Sigma*W';

P = (L'*L).*conj(R*R');
q = conj(diag(R*G'*L));
% s = trace(G'*G);

Pl = chol(P,'lower');

a = (Pl')\(Pl\q);
% a = P\q;

%%

u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end
amplitude=zeros(M,1);

for m=1:M
    aca=U*u(:,m);
    amplitude(m)=norm(aca(:),2)/sqrt(J);
end

UU=[u;deltas';omegas';amplitude']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
deltas=UU(NN+1,:);
omegas=UU(NN+2,:);
amplitude=UU(NN+3,:);
kk3=0;

for m=1:M
    if amplitude(m)/amplitude(1)>varepsilon
        kk3=kk3+1;
    else
    end
end

%% Spectral complexity: number of DMD modes.
('Spectral complexity')
kk3
u=u(:,1:kk3);
deltas=deltas(1:kk3);
omegas=omegas(1:kk3);
amplitude=amplitude(1:kk3);
('Mode number, delta, omega, amplitude')
DeltasOmegAmpl=[(1:kk3)',deltas',omegas',amplitude']

%% Reconstruction of the original snapshot matrix
hatTreconst=zeros(N,K);
for k=1:K
    hatTreconst(:,k)= ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
end

Vreconst=U*hatTreconst;

%% Calculation of DMD modes
modes=zeros(J,kk3);
amplitude0=zeros(kk3,1);
for m=1:kk3
    NormMode=norm(U*u(:,m),2)/sqrt(J);
    amplitude0(m)=NormMode;
    modes(:,m)=U*u(:,m)/NormMode;
end

%If the calculation of the amplitudes is correct, ErrAmpl=0
%ErrAmpl=norm(amplitude(:,1:kk3)-amplitude0',2)

