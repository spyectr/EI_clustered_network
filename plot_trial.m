function plot_trial(data,firings,Params)

% default plots:
plot_raster=1; % rasterplot
plot_currents=1; % postsynaptic current for representative E and I neurons
plot_stimulus=1; % stimuli time course

% load simulation results
Params.Sim.plot_length=Params.Sim.plot_length;
% save plots
savedir=Params.savedir;
if ~exist(savedir,'dir'); mkdir(savedir); end
filesave=fullfile(savedir,sprintf('SIM_[Jp%0.03g]_[p%d]',...
    Params.Jplus,Params.p));
%
Ne_plot=data.Ne_plot; % number of exc neuron to plot (typically all)
Ni_plot=data.Ni_plot; % number of inh neurons to plot  (typically all)
N_e=Ne_plot;
N_i=Ni_plot;
ind_plot=data.ind_plot;%[5; N_e+5]; % indices of neurons to plot
vplot=data.vplot; % store membrane potential for plots; rows=neurons, cols=time steps;
iEplot=data.iEplot; % store EPSC for plots; rows=neurons, cols=time steps;
iExtplot=data.iExtplot; % store IPSC for plots; rows=neurons, cols=time steps;
iIplot=data.iIplot; % store IPSC for plots; rows=neurons, cols=time steps;
VTh=data.VTh;
Jee=Params.Jee;
tau=data.tau;
p=data.p;
%
%
Ext.Mu=Params.Mu;
Network=Params.Network;
Sim=Params.Sim;
if any(strcmp(fieldnames(Params),'NcE'))
    NcE=Params.NcE;
elseif any(strcmp(fieldnames(Params),'popsize'))
    NcE=Params.popsize;
end
% keep only pops with at least 1 neuron
NcE=NcE(NcE>0);
clustermatrix=Params.clustermatrix;
clustermatrix=clustermatrix(1:numel(NcE),:);
cusumNcE=[0 cumsum(NcE)'];
fprintf('only pops sharing up to %d clusters are populated.\n',sum(clustermatrix(end,:)));

XLIM=[Sim.t_Start Sim.t_End];
if Sim.t_Start<-5
    XLIM(1)=-5;
end
if Sim.t_End>5
    XLIM(2)=5;
end

fprintf('--- Plots saved in %s...\n',filesave);
%--------------------
% FIRST PLOT: rasters
%--------------------
numfig=1;
if plot_raster
    if any(any(firings))
        figure(numfig); clf;
        numfig=numfig+1;
        % A. Plot population rasters
        option='poporder';
        fieldNames={'fieldNames','N_e','Ne_plot','Ni_plot','NcE','XLIM','option','p','clustermatrix'};
        IN=v2struct(fieldNames);
        PlotAllNeurons(firings,IN);
        % B. Ensemble rasterplot (10 units)
        saveas(gcf,[filesave '_raster_pop.pdf'],'pdf');
    end
end

%--------------------
% SECOND PLOT: PSC
%--------------------
if plot_currents
    figure(numfig); clf;
%     set(gcf, 'Position',  [100, 100, 800, 800])
    numfig=numfig+1;
    nplot=min(size(iIplot,1),2);
    np=nplot+1;
    % A. Plot membrane potentials
%         VTh=VTh*15-65;
%         vplot=vplot*15-65*ones(size(vplot)); % trasform to mV
    subplot(np,1,1);
    PlotMembrane(Sim,nplot,vplot/Jee-100,firings,ind_plot,VTh/Jee-100);
    xlim(XLIM);
    % D. Plot currents
    % bandpass IPSC
    hicutoff=200; % hi-edge (Hz)
    srate=10000; % sampling rate (Hz)
    PSCPlot(Sim,iIplot*Jee,iEplot*Jee,iExtplot*Jee,VTh*Jee./(tau),ind_plot);
    xlim(XLIM);
    % SAVE PLOT
    saveas(gcf,[filesave '_PSC.pdf'],'pdf');
end
%


%-----------------
% AUX FUN
%-----------------

% 1A. Total raster plot

function PlotAllNeurons(firings,IN)

    hold on;
    % unpack
    v2struct(IN);
    %
    cusumNcE=[0 cumsum(NcE)'];
    plot(firings(firings(:,2)<=Ne_plot,1),firings(firings(:,2)<=Ne_plot,2),'k.','markersize',1);
    plot(firings(firings(:,2)>N_e & firings(:,2)<N_e+Ni_plot,1),firings(firings(:,2)>N_e & firings(:,2)<N_e+Ni_plot,2)-N_e+Ne_plot,'r.','markersize',1);
    line([0 0],[0 Ne_plot+Ni_plot],'color','b');
    hold on;
    ylim([0 Ne_plot+Ni_plot]);
    xlab='Time (s)';
    ylab='Neurons';
    tt='All neurons ordered by population';
    figset(gca,xlab,ylab,tt,15);

    xlim(XLIM);
    hold off

end

%---------------------------------------------------

% 2. PSC Plot

function PSCPlot(Sim,iIplot,iEplot,iExtplot,VTh,ind_plot)
    stimcolors={'m',[0.5,0.5,0.5]};
    nplot=min(2,size(iIplot,1))+1;
    dt=Sim.dt_step;
    hE=[]; hI=[]; hext=[]; htot=[]; hth=[];
    t_plot=max([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
    ExcInh={'Exc','Inh'};
    Legplot={'i_E','i_I','i_{ext}','i_{tot}','V_{th}/\tau'}; hstim=[];
    for j=1:nplot-1
        subplot(nplot,1,j+1)
        hold on;
        hE(j)=plot(t_plot, iEplot(j,:),'color','b','linewidth',1); % EPSC
        hI(j)=plot(t_plot, iIplot(j,:),'color','r','linewidth',1); % IPSC
        hext(j)=plot(t_plot, iExtplot(j,:),'color','g','linewidth',1); % EPSC
        itot=iEplot(j,:)+iIplot(j,:)+iExtplot(j,:);
%         hth(j)=plot([t_plot(1) t_plot(end)],ones(1,2)*VTh(ind_plot(j))/tau(ind_plot(j)),'color','b','linewidth',1,'linestyle','-.');
        hth(j)=plot([t_plot(1) t_plot(end)],ones(1,2)*VTh(ind_plot(j)),'color','b','linewidth',1,'linestyle','-.');
        htot(j)=plot(t_plot,itot,'color','k','linewidth',1); % membrane potentials
        grid on
        xlab='Time (s)';
        ylab='mV/s';
        tt=sprintf('PSC to %s unit',ExcInh{j});
        fntsz=15;
        figset(gca,xlab,ylab,tt,fntsz);
        ind_lim=(t_plot>Sim.t_Start+0.1 & t_plot<Sim.t_End-0.1);
        legend([hE(1) hI(1) hext(1) htot(1) hth(1), hstim],Legplot);
        xlim([t_plot(1) t_plot(end)]);
        ylim([1.5*min(min(iIplot(:,ind_lim)))...
            1.5*max([max(max(iEplot(:,ind_lim))) max(max(iExtplot(:,ind_lim)))])]);
    end

end
%-------------------------------------------------

% 1A. Exc and Inh membrane potentials

function PlotMembrane(Sim,nplot,vplot,firings,ind_plot,VTh)
    ThStyle={':','--'};
    dt=Sim.dt_step;
    colors=['b'; 'r'];
    hh=[]; cnt=0;
    bins=min([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
    for j=1:nplot
        cnt=cnt+1;
        hh(cnt)=plot(bins, vplot(j,1:numel(bins)),'color',colors(j,:)); % membrane potentials
        hold on
        if any(any(firings))
            ind=find(firings(:,2)==ind_plot(j));
            if ~isempty(ind) % plot spike overshoots
                h=line([firings(ind,1), firings(ind,1)]',...
                    [VTh(ind_plot(j))*ones(size(ind)) VTh(ind_plot(j))*ones(size(ind))*2]');
                set(h,'color',colors(j,:));
                hold on
            end
        end
        cnt=cnt+1;
        hh(cnt)=plot([bins(1) bins(end)],[VTh(ind_plot(j)) VTh(ind_plot(j))],'color',colors(j,:),'linestyle',ThStyle{j});
        hold on
    end

    legend(hh,'E','V_{Th,E}','I','V_{Th,I}');
    ylim([1.2*min(min(vplot(:,0.2/dt:end))) 2*max(VTh(ind_plot))]);
    xlim([bins(1) bins(end)]);
    xlab='Time (s)';
    ylab='mV';
    tt='';
    fntsz=15;
    figset(gca,xlab,ylab,tt,fntsz);

end


%-------------------------------------------------------


function PlotPopRaster(rate,x_bins,bins)


    Cmin=min(min(rate));
    Cmax=max(max(rate));
    imagesc(x_bins,1:size(rate,1),rate); axis xy;
    caxis([Cmin, Cmax]);
    % colormap gray;
    % colormap(1-colormap);
    xlim([bins(1) bins(end)]);
    t=colorbar; get(t,'ylabel');
    set(get(t,'ylabel'),'String', 'Firing rate [spks/s]');
    hold off
    xlab='Time [s]';
    ylab='Population index';
    tt='';
    fntsz=15;
    figset(gca,xlab,ylab,tt,fntsz);

end

end