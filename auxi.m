classdef auxi
    methods(Static)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Finds Exc and Inh thresholds given firing rates in spontaneous global state
        % INPUT theta (old)
        % OUTPUT theta (new)
        %
        % from Giancarlo La Camera C_overlap_thesis
        % Luca Mazzucato February 2014
        
        
        function theta_out=fun_Fix_Threshold(ni,params)
            
            TOL_THRES=1e-9;
            TOL_POINT=1e-9;
            global DEBUG return_value screen
            
            v2struct(params);
            flag_e = 1;   % also 0: it is (almost) the same
            flag_i = 0;
            
            if DEBUG
                disp('::: threshold.c in debugging mode \n');
            end
            
            % ------------------------------------------------------------------
            %                       Starting conditions:
            % ------------------------------------------------------------------
            
            ni_e=ni(1);
            ni_i=ni(2);
            %
            ni_e_wanted = ni(1);
            ni_i_wanted = ni(2);
            
            params=auxi.fun_In_Param_2pop_v2([theta_e, theta_i],params);
            Mu=auxi.MU(ni,params); % RUN FUNCTION
            Sigma2=auxi.SIGMA2(ni,params); % RUN FUNCTION
            if return_value == 0
                fprintf('\n>>> Error in MU_C or SIGMA_C called by routine Fix_Threshold\n');
                return; % exit from function
            end
            
            %   Initial thresholds:
            theta_e = Mu(1)+3*sqrt(Sigma2(1));
            theta_i = Mu(2)+3*sqrt(Sigma2(2));
            
            %   Initial steps:
            step_e = sqrt(Sigma2(1))/10.;
            step_i = sqrt(Sigma2(2))/10.;
            
            if DEBUG
                fprintf('::: Starting Fix_Threshold...\n');
                fprintf('::: Initial parameters: Theta     [e] %5.3f  [i] %5.3f\n',theta_e,theta_i);
                fprintf('                        steps    [e] %5.3f  [i] %5.3f\n\n',step_e,step_i);
                fprintf('\n\n');
            end
            
            %------------------------------------------------------------------
            %                                 Loop:
            %------------------------------------------------------------------
            k=1; % counting loops
            while 1
                if k>20000
                    fprintf('\n>>> Error: too many steps (>%d) in routine fix_threshold \n',k-1);
                    return_value = 0;
                    return; % exit from function
                end
                
                previous_flag_e = flag_e;
                previous_flag_i = flag_i;
                params=auxi.fun_In_Param_2pop_v2([theta_e, theta_i],params);
                
                Risp=auxi.RISP(ni,params);
                if return_value == 0
                    fprintf('\n>>> Error in RISP() called by routine Fix_Threshold()\n');
                    return;
                end
                %
                ni_e = Risp(1);
                ni_i = Risp(2);
                
                if (ni_e-ni_e_wanted) > 0.0
                    flag_e = 1;
                    theta_e = theta_e+step_e;
                else
                    flag_e = 0;
                    theta_e =theta_e- step_e;
                end
                if (ni_i-ni_i_wanted) > 0.0
                    flag_i = 1;
                    theta_i = theta_i + step_i;
                else
                    flag_i = 0;
                    theta_i = theta_i- step_i;
                end
                %
                if (abs(ni_e-ni_e_wanted) < TOL_POINT && abs(ni_i-ni_i_wanted) < TOL_POINT) || (step_e< TOL_THRES && step_i< TOL_THRES)
                    if screen
                        fprintf('Threshold fixed by routine Fix_Threshold:\n');
                        fprintf('Theta_e = %g   Theta_i = %g \n\n',theta_e,theta_i);
                        fprintf('with parameters:\n');
                        fprintf('Ni_Exc: %g    Ni_Inh: %g   Ni_ext_e: %g   Ni_ext_i: %g\n\n',ni_e,ni_i,ni_ext_e,ni_ext_i);
                        fprintf('and tolerances:\n');
                        fprintf('On thresholds: %g    On fixed point: %g \n\n', TOL_THRES,TOL_POINT);
                    end
                    theta_out =[theta_e, theta_i]; % output
                    if any(strcmp(fieldnames(params),'paramsfile'))
                        save(params.paramsfile,'theta_e','theta_i','-append'); % overwrite new threholds to current paramsfile
                    end
                    return; % exit function
                end
                
                if DEBUG
                    fprintf('thresholds steps %g %g       flags %d %d      DELTA_flags %d %d  \n',step_e,step_i,flag_e,flag_i,flag_e - previous_flag_e,flag_i - previous_flag_i);
                    %fflush(stdout);
                end
                
                if (flag_e - previous_flag_e ~= 0)
                    step_e = step_e/2;          % do change!
                end
                if (flag_i - previous_flag_i ~= 0)
                    step_i = step_i/2;
                end
                
                k=k+1;
            end
            
        end
        %
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creates A and B matrices for mean and variance of the infinitesimal input
        % current to LIF (from Amit & Brunel 1997)
        function params=fun_In_Param_2pop_v2(theta,params)
            
            v2struct(params);
            
            
            % external currents
            if ~exist('ni_ext_e','var') && ~exist('ni_ext_i','var')
                ni_ext_e=ni_ext;
                ni_ext_i=ni_ext;
            end
            
            pee=Cee/N_e;
            pie=Cie/N_e;
            pei=Cei/N_i;
            pii=Cii/N_i;
            
            
            % TO AVOID PROBLEMS:
            x = 1.;
            
            A(1,1) = tau_e*Cee*x*Jee;
            A(1,2) = -tau_e*Cei*Jei;
            A(2,1) = tau_i*Cie*x*Jie;
            A(2,2) = -tau_i*Cii*Jii;
%             
%             % variance for fixed in degree C (like in Amit-Brunel 1997)
%             B(1,1) = tau_e*Cee*x*Jee*Jee;
%             B(1,2) = tau_e*Cei*Jei*Jei;
%             B(2,1) = tau_i*Cie*x*Jie*Jie;
%             B(2,2) = tau_i*Cii*Jii*Jii;
            
            % variance for Erdos-Renyi (like in Renart et al Science 2010)
            B(1,1) = tau_e*N_e*pee*(1-pee)*x*Jee*Jee;
            B(1,2) = tau_e*N_e*pei*(1-pei)*Jei*Jei;
            B(2,1) = tau_i*N_i*pie*(1-pie)*x*Jie*Jie;
            B(2,2) = tau_i*N_i*pii*(1-pii)*Jii*Jii;
            
            
            
            if (delta ~= 0)
                for i=1:2
                    for j=1:2 B(i,j) = B(i,j) *(1.+delta^2); end
                end
            end
            
            % if (x!=1)  Cext = C*(1-x);
            
            Mu_ext(1) = Cext*ni_ext_e*tau_e*Jee_ext;
            Mu_ext(2) = Cext*ni_ext_i*tau_i*Jie_ext;
            
            Sigma_ext=zeros(1,2);
%             if strcmp(Stimulus.input,'Poisson')
%                 Sigma_ext(1) = (1.+delta^2)*Cext*ni_ext_e*tau_e*Jee_ext^2;
%                 Sigma_ext(2) = (1.+delta^2)*Cext*ni_ext_i*tau_i*Jie_ext^2;
%             end
            
            H(1) = He*theta(1);
            H(2) = Hi*theta(2);
            
            Tau(1) = tau_e;
            Tau(2) = tau_i;
            
            Tausyn(1) = tausyn_e;
            Tausyn(2) = tausyn_i;
            
            % params=[]; % output structure
            params.A=A;
            params.B=B;
            params.Mu_ext=Mu_ext;
            params.Sigma_ext=Sigma_ext;
            params.H=H' ;
            params.Tau=Tau';
            params.Tausyn=Tausyn';
            params.Theta=theta';
            params.tau_arp=tau_arp;
%             % vector d (dimensions) is defined in _espo only if overlap>0; for the
%             % moment let's define it as a zeros:
%             params.d=zeros(1,nn);
            % synaptic dynamics
%             if strcmp(Network.syn,'DoubleExp') || strcmp(Network.syn,'SingleExp')
                params.BS=1; % activates Brunel-Sergi thresholds
%             else
%                 params.BS=0; % usual transfer function
%             end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transfer function
        % INPUT ni   = firing rate row vector
        %       Params = input parameters
        %       Theta (optional): if absent, thresholds are set from Params.Theta
        %                         if present, thresholds are set equal to Theta
        % OUTPUT Risp=firing rate vector;
        
        %
        % from Giancarlo La Camera C_overlap_thesis
        % Luca Mazzucato February 2014
        % Luca Mazzucato May 2020
        
        function Risp=RISP(ni,Params)
            
            warning('off','all');
            
            global return_value
            
            n=numel(ni); % # of active populations
            
            % QUENCH=1: option for quenched noise
            QUENCH=0;
            if any(strcmp(fieldnames(Params),'Network'))
                if any(strcmp(fieldnames(Params.Network),'cuemode'))
                    %         if strcmp(Params.Network.cuemode,'gaussian') && n>2
                    if ~isempty(strfind(Params.Network.cuemode,'gaussian')) && n>2
                        % if n=2 override gaussian cue (only used for threshold
                        % computations)
                        QUENCH=1;
                        quench_pops=Params.quench_pops;
                        nonquench_pops=setxor(1:n,quench_pops);
                    end
                end
            end
            
            % parameters
            H=Params.H;
            Tau=Params.Tau;
            tau_arp=Params.tau_arp;
            Theta=Params.Theta;
            
            Risp=zeros(n,1);
            % k = 0.9;
            % Brunel Sergi and Brunel Fourcaud: correction to transfer function to
            % account for exponential synapses
            Tausyn=Params.Tausyn;
            BS=zeros(n,1);
            if Params.BS==1
                a=1.032626576115609; %-zeta(1/2)/sqrt(2); %
                BS=a*sqrt(Tausyn./Tau);
            end
            
            Mu=auxi.MU(ni,Params);         % calcola le Mu(ni)
            
            Sigma=auxi.SIGMA2(ni,Params);      % calcola le Sigma quadro(ni)
            if(return_value == 0)
                fprintf('\n>>> Error in routine SIGMA() called by RISP()\n');
                return;
            end
            if any(Sigma<0)
                fprintf('\n SIGMA<0!!\n ');
            end
            if QUENCH==0
                % usual transfer function
                for i=1:n
                    sigma = sqrt(Sigma(i));
                    thres = (Theta(i)-Mu(i))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset = (H(i)-Mu(i))/sigma+BS(i);
                    % NOTE: using Qsimp or quad, the integral gives the wrong values using nerf
                    % for w<0, so we need to use Maurizio Mattia's trick to get the correct
                    % numerical result. However, it takes 20x w.r.t. the mex file.
                    %     integ = Qsimp(@auxPhiExp, reset, thres);
                    %     integ = quad(@auxPhiExp, reset, thres);
                    %     integ=1.772453851*integ;
                    % USE MEX FILE TO EVALUATE INTEGRAL (already includes sqrt(pi) factor)
                    % NOTE: the MEX file integrates nerf both for w>=0 and w<0, giving in
                    % both cases the correct result without Mattia's trick!
                    integ=IntAuxPhi_vec(reset,thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine qsimp() called by RISP()\n');
                        return_value = 0;
                        return;
                    end
                    
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp(i) = 1./inv_risp;
                end
            elseif QUENCH==1
                % quenched gaussian noise on clustered populations only
                for i=quench_pops
                    sigma = sqrt(Sigma(i));
                    MuZ=@(z)(Mu(i)+Params.Mu_extZ(z)); % input current with gaussian quenched noise z
                    thres =@(z) (Theta(i)-MuZ(z))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset =@(z) (H(i)-MuZ(z))/sigma+BS(i);
                    integfun=@(z)(exp(-z.^2/2)/sqrt(2*pi))./(tau_arp + Tau(i)*IntAuxPhi_vec(reset(z),thres(z)));
                    integ=integral(integfun,-20,20);
                    if (integ < 0)
                        fprintf('\n>>> Error in integral of IntAuxPhi_vec called by RISP_DYN()\n');
                        return_value = 0;
                        return;
                    end
                    Risp(i)=integ;
                end
                % background and inhibitory populations without quenched noise
                for i=nonquench_pops
                    sigma = sqrt(Sigma(i));
                    thres = (Theta(i)-Mu(i))/sigma+BS(i);
                    %thres = (Theta(i)+(k*ni(i))-Mu(i))/sigma;% Calcium effect
                    reset = (H(i)-Mu(i))/sigma+BS(i);
                    integ=IntAuxPhi_vec(reset,thres);
                    if (integ == -1)
                        fprintf('\n>>> Error in routine qsimp() called by RISP()\n');
                        return_value = 0;
                        return;
                    end
                    inv_risp = (tau_arp + Tau(i)*integ);
                    Risp(i) = 1./inv_risp;
                end
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        
        function Mu=MU(ni,Params)
            
            % global return_value
            
            nn=numel(ni);
            A=Params.A;
            Mu_ext=Params.Mu_ext;
            Mu_perc = zeros(nn+1);
            Mu=zeros(nn,1);
            
            
            for i=1:nn
                Mu(i) = Mu_ext(i);
                for j=1:nn
                    Mu(i)=Mu(i)+ A(i,j)*ni(j);
                end
            end
            
            for i=1:nn
                for j=1:nn
                    % percentage of afferent current j->i
                    Mu_perc(i,j) = A(i,j)*ni(j);
                end
                % percentage of afferent external current
                Mu_perc(i,nn+1) = Mu_ext(i);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % from GLC C_overlap_thesis
        % Luca Mazzucato February 2014
        
        
        function Sigma2=SIGMA2(ni,Params)
            
            global return_value
            
            B=Params.B;
            Sigma_ext=Params.Sigma_ext;
            
            Sigma2=B*ni'+Sigma_ext';
            
            if any(Sigma2<0.0)
                
                fprintf('\n\n>>> Error: Sigma[%d]<0\n',find(Sigma2<0));
                return_value = 0;
                return;
            end
        end
        
        
    end
end