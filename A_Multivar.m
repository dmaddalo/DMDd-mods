clearvars; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE NUMBER 3 / MULTI VARIABLE
% This code takes the variables saved previously and performs the
% algorithms on the selected variable TOGETHER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% VARIABLES
% ELECTRIC POTENTIAL = phi
% PLASMA DENSITY = n
% ELECTRON TEMPERATURE = Te
% ION AXIAL VELOCITY = vi1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set as true only if you still have to figure out best tolerances values,
% otherwise keep this flag zero all the time
tolerances = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc = true;
ca = 3;
va = {'n','Te','phi','vi1'};                        % KEEP THIS ORDER HERE
if tolerances == false
    eps1 = 5e-3;
    eps11 = 5e-2;
    eps2 = 3e-2;
end
d = 700;
skip = 1;
cutt = 0;

PLOTSVD = 0;
PLOTDMD = 0;
SAVESVD = 0;
SAVEDMD = 0;
RECODMD = 0;                                                                                                                                                                                                                                                                                                                                                                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full matrix construction.
for j = 1:size(va,2)
    if acc == false
        load(['sim',num2str(ca),'_',char(va(j)),'3d_be.mat'],'n3d');
    else
        load(['sim',num2str(ca),'_',char(va(j)),'3d_be_acc.mat'],'n3d');
    end
    
    % n3d = flip(n3d,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Snapshot cut-off for possible unwanted transients.
    n3d = n3d(:,:,cutt+1:end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Carrying all snapshots is not always necessary. Each
    % dataset/variable may require a different "skip" value.
    n3d = n3d(:,:,1:skip:end);
    
    mn3d = mean(n3d,3); fn3d = n3d - mn3d;
    
    Vp = zeros(size(n3d,1)*size(n3d,2),size(n3d,3));
    for i = 1:size(n3d,3)
        Vp(:,i) = reshape(fn3d(:,:,i),[1 size(n3d,1)*size(n3d,2)]);
    end
    ss = sqrt(sum(vecnorm(Vp,2).^2)); Vp = Vp/ss;
    
    if j == 1
        V = Vp;
    else
        V = [V;Vp];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same cut-off and skipping procedure for time vector.
load(['sim',num2str(ca),'_t_be.mat'],'t');
t = t(cutt+1:end);
t = t(1:skip:end);

[My,Mx,N] = size(n3d);
dt = t(2) - t(1); fsam = 1/dt; fnyq = fsam/2;
fRFT = linspace(0,fnyq,(N-1)/2+1);

clear Vp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVD RESULTS
% Select and plot the number of desired POD modes with relative time
% envelope and Fourier spectrum. Each window contains one mode. Turn
% "PLOTSVD" on above in order to see it.
mds = 11;

cmap = [0,0,0; jet];

[U,Sfull,W] = svd(V,'econ');

if PLOTSVD == true
    for j = 1:mds
        figure
        tiledlayout('flow')
        for i = 1:4
            nexttile
            mode = reshape(U((Mx*My)*(i-1)+1:(Mx*My)*i,j),[My Mx]);
            modev = std(mode,0,'all'); modem = mean(mode,'all');
            mode(29:41,1:19) = NaN; mode(1:11,1:19) = NaN;
            imagesc(mode); colormap(cmap); colorbar;
            caxis([modem-4*modev modem+4*modev]);
            axis square;
            switch i
                case 1
                    title('plasma density');
                case 2
                    title('electron temperature');
                case 3
                    title('plasma potential');
                case 4
                    title('ion axial velocity');
            end
        end

        FT = abs(fft(W(:,j)));
        RFT = FT(1:(size(FT,1)+1)/2);
        nexttile
        plot(t,W(:,j))
        title('time domain');
        nexttile
        semilogx(fRFT,RFT); axis([-inf fnyq -inf inf]);
        title('fourier domain'); set(gca,'Xgrid','on');

        sgtitle(['mode ',num2str(j)]);

        set(gcf,'Units','Normalized','OuterPosition',[-0.001,0.035,1.002,0.96]);
    end
end

clear modem modev

if SAVESVD == true
    ifinal = get(gcf,'Number');
    for i = 1:ifinal
        i
        if acc == false
            saveas(figure(i),[pwd '/AdriansFigures/sim',num2str(ca),'_svdmode',num2str(i),'_inst.png']);
        else
            saveas(figure(i),[pwd '/AdriansFigures/sim',num2str(ca),'_svdmode',num2str(i),'_acc.png']);
        end
    end
end

% SVD SPECTRUM

% figure
% line(1:size(Sfull,2),diag(Sfull)/Sfull(1,1));
% set(gca,'xscale','log','yscale','log');
% xlabel('mode'); ylabel('\sigma_i/\sigma_1');
% title('SVD spectrum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DMD RESULTS
% Keep "eps1" and "eps11" the same. Each variable may require different
% values of tolerance due to subjective error levels. Multivariable case
% requires an energy cut-off of 1e-2, not below. Spectra are plot both in
% logarithmic scale for both axes and symmetrically with respect to x axis.

if tolerances == true
    [eps1,eps11] = A_Tol(d,V);
end

[Vrecon,deltas,omegas,amplitudes,modes] = DMDd_SIADS(d,V,t,eps1,eps11,eps2);

s_deltas = sign(deltas);
unst = abs(deltas*dt) > 0.69/N;

figure;
tiledlayout('flow')
%% XSYMMETRIC PLOT 1
nexttile;
scatter(omegas(s_deltas == 1)/(2*pi),deltas(s_deltas == 1)*dt); hold on;
scatter(omegas(s_deltas == -1)/(2*pi),abs(deltas(s_deltas == -1))*dt);
xlabel('Frequency (Hz)'); ylabel('Growth rates'); yline(0.69/N);
scatter(omegas(unst)/(2*pi),abs(deltas(unst))*dt,'*');
axis([-2e5 2e5 1e-6 1e0]);
set(gca,'yscale','log');
grid on; grid minor;
legend('Positive growth rates','Negative growth rates',...
    'Threshold (double/half)','location','best');
%% XSYMMETRIC PLOT 2
nexttile;
scatter(omegas(s_deltas == 1)/(2*pi),amplitudes(s_deltas == 1)); hold on;
scatter(omegas(s_deltas == -1)/(2*pi),amplitudes(s_deltas == -1));
xlabel('Frequency (Hz)'); ylabel('Amplitudes');
scatter(omegas(unst)/(2*pi),amplitudes(unst),'*');
axis([-2e5 2e5 1e-5 5e-3]);
set(gca,'yscale','log');
grid on; grid minor;
legend('Positive growth rates','Negative growth rates','location','best');
%% LOGLOG PLOT 1
nexttile;
scatter(omegas(s_deltas == 1)/(2*pi),deltas(s_deltas == 1)*dt); hold on;
scatter(omegas(s_deltas == -1)/(2*pi),abs(deltas(s_deltas == -1))*dt);
xlabel('Frequency (Hz)'); ylabel('Growth rates'); yline(0.69/N);
scatter(omegas(unst)/(2*pi),abs(deltas(unst))*dt,'*');
axis([1e3 fnyq 1e-6 1e0]);
set(gca,'yscale','log','xscale','log');
grid on; grid minor;
legend('Positive growth rates','Negative growth rates',...
    'Threshold (double/half)','location','best');
%% LOGLOG PLOT 2
nexttile;
scatter(omegas(s_deltas == 1)/(2*pi),amplitudes(s_deltas == 1)); hold on;
scatter(omegas(s_deltas == -1)/(2*pi),amplitudes(s_deltas == -1));
xlabel('Frequency (Hz)'); ylabel('Amplitudes');
scatter(omegas(unst)/(2*pi),amplitudes(unst),'*');
axis([1e3 fnyq 1e-5 5e-3]);
set(gca,'yscale','log','xscale','log');
grid on; grid minor;
legend('Positive growth rates','Negative growth rates','location','best');

if acc == 0
    sgtitle({['d = ',num2str(d),'; skip = ',num2str(skip),'; cut = ',num2str(cutt),...
    '; inst variables'],['Energy tolerance = ',sprintf('%0.1e',eps1)],...
    ['Amplitude tolerance = ',sprintf('%0.1e',eps2*max(amplitudes)),...
    ' (',num2str(eps2*100),'%)']});
else
    sgtitle({['d = ',num2str(d),'; skip = ',num2str(skip),'; cut = ',num2str(cutt),...
        '; acc variables'],['Energy tolerance = ',sprintf('%0.1e',eps1)],...
        ['Amplitude tolerance = ',sprintf('%0.1e',eps2*max(amplitudes)),...
        ' (',num2str(eps2*100),'%)']});
end
set(gcf,'Position',[10 40 1120 840]);
% set(gcf,'Position',[10 40 560 840]);
%% END OF HODMD SPECTRA PLOTS

modeamps = abs(modes); modephases = angle(modes);

if PLOTDMD == true
        for i = 1:2:size(modes,2)
            figure
            p = uipanel('Title',['Freq: ',sprintf('%0.2e',omegas(i)/(2*pi)),...
                '; GR: ',sprintf('%0.2e',deltas(i)*dt),'; Amplitude',sprintf('%0.2e',amplitudes(i)),...
                '; Modulus'],'Position',[0 0.5 1 0.5]);
            tiled = tiledlayout(p,1,4);
            for j = 1:4
                nexttile(tiled)
                mode = reshape(modeamps((Mx*My)*(j-1)+1:(Mx*My)*j,i),[My Mx]);
                modev = std(mode,0,'all'); modem = mean(mode,'all');
                mode(29:41,1:19) = NaN; mode(1:11,1:19) = NaN;
                imagesc(mode); colormap(cmap); colorbar;
                caxis([modem-4*modev modem+4*modev]);
                axis square;
                switch j
                    case 1
                        title('Plasma density');
                    case 2
                        title('Electron temperature');
                    case 3
                        title('Plasma potential');
                    case 4
                        title('Ion axial velocity');
                end
            end
            p = uipanel('Title',['Freq: ',sprintf('%0.2e',omegas(i)/(2*pi)),...
                '; GR: ',sprintf('%0.2e',deltas(i)*dt),'; Amplitude',sprintf('%0.2e',amplitudes(i)),...
                '; Phase'],'Position',[0 0 1 0.5]);
            tiled = tiledlayout(p,1,4);
            for j = 1:4
                nexttile(tiled)
                mode = reshape(modephases((Mx*My)*(j-1)+1:(Mx*My)*j,i),[My Mx]);
                modev = std(mode,0,'all'); modem = mean(mode,'all');
                mode(29:41,1:19) = 0; mode(1:11,1:19) = 0;
                imagesc(mode); set(gca,'colormap', hsv); colorbar;
                caxis([-pi pi]);
                axis square;
                switch j
                    case 1
                        title('Plasma density');
                    case 2
                        title('Electron temperature');
                    case 3
                        title('Plasma potential');
                    case 4
                        title('Ion axial velocity');
                end
            end
            set(gcf,'Units','Normalized','OuterPosition',[-0.001,0.035,1.002,0.96]);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DMD RECONSTRUCTION
% The reconstruction is provided into a cell array "Recons", which number
% of cells corresponds to the number of clusters selected. "Recons" is then
% redistributed into different cell arrays, one for each variable, having
% the second row as the standard deviation used for having a good scaling
% on the colorbar. Turn on the "RECONDMD" flag to see the reconstruction.

if RECODMD == 1
    Recons = A_Recon(deltas,omegas,amplitudes,modes,t);

    clust = size(Recons,2);

    n3dR = cell(2,clust); Te3dR = n3dR; phi3dR = n3dR; vi13dR = n3dR;
    for j = 1:clust
        for i = 1:N
            n3dR{1,j}(:,:,i) = reshape(Recons{1,j}(1:My*Mx,i),[My,Mx,1]);
            n3dR{2,j}(i) = std(Recons{1,j}(1:My*Mx,i),0,1);
            Te3dR{1,j}(:,:,i) = reshape(Recons{1,j}(My*Mx+1:My*Mx*2,i),[My,Mx,1]);
            Te3dR{2,j}(i) = std(Recons{1,j}(My*Mx+1:My*Mx*2,i),0,1);
            phi3dR{1,j}(:,:,i) = reshape(Recons{1,j}(2*My*Mx+1:My*Mx*3,i),[My,Mx,1]);
            phi3dR{2,j}(i) = std(Recons{1,j}(2*My*Mx+1:My*Mx*3,i),0,1);
            vi13dR{1,j}(:,:,i) = reshape(Recons{1,j}(3*My*Mx+1:My*Mx*4,i),[My,Mx,1]);
            vi13dR{2,j}(i) = std(Recons{1,j}(3*My*Mx+1:My*Mx*4,i),0,1);
        end
        n3dR{2,j} = mean(n3dR{2,j}); Te3dR{2,j} = mean(Te3dR{2,j});
        phi3dR{2,j} = mean(phi3dR{2,j}); vi13dR{2,j} = mean(vi13dR{2,j});
    end

    cl = 2; see =vi13dR;

    figure
    for i = 1:N
        see{1,cl}(29:41,1:19,i) = NaN; see{1,cl}(1:11,1:19,i) = NaN;
        imagesc(see{1,cl}(:,:,i)); 
        colormap(cmap); colorbar; text(40,4,num2str(t(i)),'Color','w')
        caxis([-10*see{2,cl} 10*see{2,cl}]); pause(0.01);
    end
end

