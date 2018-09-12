function [fitMM_ls_mle,res]=...
    fitMM_ls_mle(sig,scanParam,startingVals,lwrBnds,...
    upprBnds,plot_y_n,chooseModel,ls_mle_str,noise)
%% Function for fitting microstructural model with either:
%  one diffusivitiy, D; with or without an fiso compartment
%  two diffusivities, Di & De; with or without an fiso compartment
%% Inputs: sig - vector of signals
%          scanParam - matrix of scan parameters
%          startingValus - matrix of starting values to use in fitting
%          lwrBnd/upprBnd - vector of lower and upper bounds:
%                         use  lwrBnd=[-Inf -Inf -Inf -Inf] and
%                              upprBnd=[Inf Inf Inf Inf] for unconstrained
%                              fitting
%                         use  lwrBnd=[-Inf x -Inf -Inf] and
%                              upprBnd=[Inf x Inf Inf] to fix second
%                              parameter at x
%        plot_y_n - string indicating whether to plot ('y') or not ('n')
%        chooseModel - either 'D', 'D_fiso', 'DiDe', 'DiDe_fiso'
%        ls_mle_str - 'mle' to specify non-ls fitting obj fn
%        noise - Rician noise SD vector

%% Outputs: fitMM_ls_mle - fitted parameter estimates
%           res - residuals from fit

%% Load sphere GPA equation roots - Neuman, J Chem Phys 1974;60:4508-4511.
rootsS=...
	loadRename('./+mm_fitting/gpaBesselEqnRoots_sphere.mat');

%% Fit model and collect parameter estimates
options=optimset('Display', 'off');

%% Check if we're fitting with LS or MLE
switch ls_mle_str
    case 'mle'
        ls_mle='mle';
        riceNoiseStd=noise;
        %disp('*********************************')
        %disp('******TESTING MLE FITTING******')
        %disp('*********************************')
    otherwise
        ls_mle='ls';
end

%% Loop over initial values and accept fit which gives lowest value of objective function
ofv(1)=Inf;
for strtValInd=1:size(startingVals,2)
    %{%
    %SIMPLEX:
    switch chooseModel
        case {'D','DiDe','DiDe_notort','DiDe_with_rise_time'}
            [fitParams,objFnVal,extflg,output]=mm_fitting.fminsearchbnd(@subFn,...
                startingVals(:,strtValInd),lwrBnds,upprBnds,options);
        case 'DiDe_unnormalised'
            % include starting values for S0
            g0MeanSig=mean(sig(find(scanParam(:,1)==0)));
            startingValsToUse=startingVals(:,strtValInd);
            startingValsToUse(end)=g0MeanSig;
            [fitParams,objFnVal,extflg,output]=mm_fitting.fminsearchbnd(@subFn,...
                startingValsToUse,lwrBnds,upprBnds,options);
        case 'DiDe_fiso'  % constrain fi+fe<1
            [fitParams,objFnVal,extflg,output]=mm_fitting.fminsearchcon(@subFn,...
                startingVals(:,strtValInd),lwrBnds,upprBnds,[0 0 0 1 1],[1],[],options);
            %options=optimoptions(@fmincon,'Algorithm','sqp',...
            %    'Display','off');
            %[fitParams,objFnVal,extflg,output]=fmincon(@subFn,...
            %    startingVals(:,strtValInd),[0 0 0 1 1],1,[],[],...
            %    lwrBnds,upprBnds,[],options);
            %[fitParams,objFnVal,extflg]=nag_e04wd_wrapper(@subFnNag_e04wd,...
            %    startingVals(:,strtValInd),[0 0 0 1 1],[0],[1],...
            %    lwrBnds,upprBnds,[],'cold');
            %[fitParams,objFnVal,extflg]=nag_e04us_wrapper(@subFnNag_e04us,...
            %    sig,startingVals(:,strtValInd),[0 0 0 1 1],[0],[1],...
            %    lwrBnds,upprBnds,[],'cold');
        case 'DiDe_fit_fiso'  % fmincon needed to constrain fi+fe+fiso=1
            options=optimoptions(@fmincon,'Algorithm','sqp',...
                'Display','off');
            options=optimoptions(@fmincon,'Display','off');
            [fitParams,objFnVal,extflg,output]=fmincon(@subFn,...
                startingVals(:,strtValInd),[],[],[0 0 0 1 1 1],1,...
                lwrBnds,upprBnds,[],options);
            %[fitParams,objFnVal,extflg]=nag_e04wd_wrapper(@subFnNag,...
            %    startingVals(:,strtValInd),[0 0 0 1 1 1],[1],[1],...
            %    lwrBnds,upprBnds,[],'cold');
            %[fitParams,objFnVal,extflg]=nag_e04us_wrapper(@subFnNag_e04us,...
            %    sig,startingVals(:,strtValInd),[0 0 0 1 1 1],[1],[1],...
            %    lwrBnds,upprBnds,[],'cold');
            % Try genetic algorithm
            %options = gaoptimset('InitialPopulation',startingVals',...
            %    'MutationFcn',@mutationadaptfeasible,'Display','off',...
            %    'StallGenLimit',100,'PopulationSize',size(startingVals,2));
            %'PlotFcns',@gaplotbestf,'Vectorized','on');
            %[fitParams,objFnVal,extflg,output,population]=...
            %    ga(@subFn,6,[],[],...
            %    [0 0 0 1 1 1],[1],lwrBnds,upprBnds,[],options);
            %
        case 'D_fiso' % constrain fi+fe<1
            [fitParams,objFnVal,extflg,output]=mm_fitting.fminsearchcon(@subFn,...
                startingVals(:,strtValInd),lwrBnds,upprBnds,[0 0 1 1 0],[1],[],options);
    end
    
    ofv(strtValInd+1)=objFnVal;
    % Update estimates if obj fun is lower
    if ofv(strtValInd+1)<min(ofv(1:strtValInd))
        switch chooseModel
            case 'D'
                rFit=fitParams(1);
                dFit=fitParams(2);
                fFit=fitParams(3);
                s0Fit=fitParams(4);
                
            case 'D_fiso'
                rFit=fitParams(1);
                dFit=fitParams(2);
                fiFit=fitParams(3);
                feFit=fitParams(4);
                s0Fit=fitParams(5);
                extFlg=extflg;
                
            case {'DiDe','DiDe_with_rise_time'}
                rFit=fitParams(1);
                diFit=fitParams(2);
                deFit=fitParams(3);
                fFit=fitParams(4);
                
            case 'DiDe_unnormalised'
                rFit=fitParams(1);
                diFit=fitParams(2);
                deFit=fitParams(3);
                fFit=fitParams(4);
                s0Fit=fitParams(5);
                
            case 'DiDe_notort'
                rFit=fitParams(1);
                diFit=fitParams(2);
                deFit=fitParams(3);
                fFit=fitParams(4);
                
            case 'DiDe_fiso'
                rFit=fitParams(1);
                diFit=fitParams(2);
                deFit=fitParams(3);
                fiFit=fitParams(4);
                feFit=fitParams(5);
                extFlg=extflg;
                
            case 'DiDe_fit_fiso'
                rFit=fitParams(1);
                diFit=fitParams(2);
                deFit=fitParams(3);
                fiFit=fitParams(4);
                feFit=fitParams(5);
                fisoFit=fitParams(6);
                extFlg=extflg;
                
            otherwise
                error('!!! Invalid option for chooseModel')
        end
        strtValIndAccepted=strtValInd;
    else
        % If these starting values do not improve obj fn, do nothing
    end
end

switch chooseModel
    case 'D'
        fitParamsAccepted=cat(2,rFit,dFit,fFit,s0Fit);
    case 'D_fiso'
        fitParamsAccepted=cat(2,rFit,dFit,fiFit,feFit,s0Fit);
    case {'DiDe', 'DiDe_notort','DiDe_with_rise_time'}
        fitParamsAccepted=cat(2,rFit,diFit,deFit,fFit);
    case 'DiDe_unnormalised'
        fitParamsAccepted=cat(2,rFit,diFit,deFit,fFit,s0Fit);    
    case 'DiDe_fiso'
        fitParamsAccepted=cat(2,rFit,diFit,deFit,fiFit,feFit);
    case 'DiDe_fit_fiso'
        fitParamsAccepted=cat(2,rFit,diFit,deFit,fiFit,feFit,fisoFit);
end

%% Calculate R^2 of fit
%  calculate residuals for each data point
res = sig - subFnPlot(fitParamsAccepted,scanParam);
sumSqRes = sum(res.^2);
diffFromMean = sig - mean(sig);
sumSqDiffFromMean = sum(diffFromMean.^2);
rsq = 1 - (sumSqRes./sumSqDiffFromMean);

%% Assign output
switch chooseModel
    case {'DiDe_fit_fiso','DiDe_fiso','D_fiso'}
        fitMM_ls_mle=...
            cat(1,fitParamsAccepted(:),rsq,strtValIndAccepted,extFlg);
    case {'D','DiDe','DiDe_notort','DiDe_unnormalised','DiDe_with_rise_time'}
        fitMM_ls_mle=...
            cat(1,fitParamsAccepted(:),rsq,strtValIndAccepted);
end

%% Plotting, if required
switch plot_y_n
    case 'y'
        DELvals=unique(scanParam(:,2));
        figure, hold on
        C= [0    0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];
        for i=1:size(DELvals(:),1)
            indexToDEL=DELvals(i);%sequenceParams(indexToDEL,2)
            [r,~,~]=find(scanParam(:,2)==DELvals(i));
            G=scanParam(r,1);
            xgrid = linspace(min(G),max(G))';
            plot(G.*1e3,sig(r),'o',xgrid.*1e3,subFnPlot(fitParamsAccepted,cat(2,xgrid,...
                repmat(scanParam(r(1),2:end),size(xgrid,1),1))),'-.',...
                'MarkerSize',10,'LineWidth',2,'Color',C(i,:));
        end
        xlabel('G (mT/m)','FontSize',22);
        ylabel('Signal intensity','FontSize',22);
        ylim([0 1])
        set(gca,'YGrid','on','box','off','FontSize',22);
        
        % Add text displaying fitted values
        switch chooseModel
            case 'D'
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nD=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(3)),...
                    '\nS0=',num2str(fitMM_ls_mle(4)),...
                    '\nR^2=',num2str(fitMM_ls_mle(5)))),...
                    'fontsize',16);
            case 'D_fiso'
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nD=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(3)),...
                    '\nfe=',num2str(fitMM_ls_mle(4)),...
                    '\nfiso=',num2str(1-(fitMM_ls_mle(4) + ...
                    fitMM_ls_mle(5))),...
                    '\nS0=',num2str(fitMM_ls_mle(5)),...
                    '\nR^2=',num2str(fitMM_ls_mle(6)))),...
                    'fontsize',16);
            case {'DiDe','DiDe_notort','DiDe_with_rise_time'}
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nDi=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nDe=',num2str(fitMM_ls_mle(3)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(4)),...
                    '\nR^2=',num2str(fitMM_ls_mle(5)))),...
                    'fontsize',16);
            case 'DiDe_unnormalised'
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nDi=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nDe=',num2str(fitMM_ls_mle(3)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(4)),...
                    '\ns0=',num2str(fitMM_ls_mle(5)),...
                    '\nR^2=',num2str(fitMM_ls_mle(6)))),...
                    'fontsize',16);    
            case 'DiDe_fiso'
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nDi=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nDe=',num2str(fitMM_ls_mle(3)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(4)),...
                    '\nfe=',num2str(fitMM_ls_mle(5)),...
                    '\nfiso=',num2str(1-(fitMM_ls_mle(4) + ...
                    fitMM_ls_mle(5))),...
                    '\nR^2=',num2str(fitMM_ls_mle(6)))),...
                    'fontsize',16);
            case 'DiDe_fit_fiso'
                text(0.15*max(G.*1e3),0.3*max(sig),sprintf(strcat(...
                    'R=',num2str(fitMM_ls_mle(1)*1e6),' \\mum',...
                    '\nDi=',num2str(fitMM_ls_mle(2)*1e9), '\\mum^2/ms',...
                    '\nDe=',num2str(fitMM_ls_mle(3)*1e9), '\\mum^2/ms',...
                    '\nfi=',num2str(fitMM_ls_mle(4)),...
                    '\nfe=',num2str(fitMM_ls_mle(5)),...
                    '\nfiso=',num2str(fitMM_ls_mle(6)),...
                    '\nR^2=',num2str(fitMM_ls_mle(7)))),...
                    'fontsize',16);
        end
    case 'n'
        %do nothing
    otherwise
        error('!!! plot_y_n must be either y or n')
end

%% Signal model functions
% Sphere expression - Murday & Cotts, J Chem Phys 1968;48:4938–4945.
% Extracellular tortuosity - Price et al., Biophys J 1998;74:2259–2271.
    function x=subFn(a)
        %sumM=0;
        Grad=scanParam(:,1);
        DEL=scanParam(:,2);
        del=scanParam(:,3);
        gamma=scanParam(:,4);
        switch chooseModel
            case 'D'
                R=a(1);
                D=a(2);
                f=a(3);
                S0=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                
            case 'D_fiso'
                R=a(1);
                D=a(2);
                fi=a(3);
                fe=a(4);
                S0=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                
            case {'DiDe','DiDe_unnormalised'}
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                switch chooseModel
                    case 'DiDe'
                        % nothing to do
                    case 'DiDe_unnormalised'
                        s0=a(5);
                end
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                switch chooseModel
                    case 'DiDe'
                        signalModel=...
                            1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                    case 'DiDe_unnormalised'
                        signalModel=...
                            s0.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                end
            
            case 'DiDe_with_rise_time'
                %
                e=scanParam(:,5);
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                sumTerm = 0;
                N = 1; % Number of half oscillations (i.e. N = 1 for PGSE)
                v = N./(2.*del);
                for n=1:numel(rootsS)
                    u = rootsS(n);
                    l = (u/R)^2;
                    B = (2*(R/u)^2)/(u^2 - 2);
                    % Waveform
                    % Gamma2 in GPD_trapezoid_arbitrary_freq.m - using ToMatlab package
                    W=(-2).*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^(-1).*(N+2.*(del+DEL).* ...
                        v)).*(exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)) ...
                        ).^4.*Grad.^2.*l.^(-2)+2.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^( ...
                        -1).*(N+2.*(2.*e+DEL).*v)).*((-1)+exp(1).^(Di.*e.*l)).^2.*(exp(1).^(Di.* ...
                        e.*l)+(-1).*exp(1).^((1/2).*Di.*l.*v.^(-1))).^2.*(1+exp(1).^((1/2).* ...
                        Di.*l.*v.^(-1))).^(-2).*((-1).^(1+N)+exp(1).^((1/2).*Di.*N.*l.*v.^(-1))) ...
                        .*((-1)+(-1).^N.*exp(1).^((1/2).*Di.*N.*l.*v.^(-1))).*Grad.^2.*l.^(-2)+2.*( ...
                        -1).^N.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^(-1).*(N+2.*(e+ ...
                        del+DEL).*v)).*((-1)+exp(1).^(Di.*e.*l)).*(exp(1).^(Di.*e.*l)+(-1).*exp(1) ...
                        .^((1/2).*Di.*l.*v.^(-1))).*(1+exp(1).^((1/2).*Di.*l.*v.^(-1))).^(-1).*( ...
                        exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1))).^2.* ...
                        (1+(-1).^(1+N).*exp(1).^(Di.*del.*l)+(-2).*exp(1).^(Di.*DEL.*l)+exp(1).^( ...
                        Di.*del.*l+(1/2).*Di.*N.*l.*v.^(-1))+2.*(-1).^N.*exp(1).^(Di.*DEL.*l+(1/2).* ...
                        Di.*N.*l.*v.^(-1))+(-1).^(1+N).*exp(1).^((1/2).*Di.*N.*l.*v.^(-1))).* ...
                        Grad.^2.*l.^(-2)+(-4).*(-1).^N.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.* ...
                        l.*v.^(-1).*(1+N+4.*e.*v)).*((-1)+exp(1).^(Di.*e.*l)).^2.*(exp(1).^( ...
                        Di.*e.*l)+(-1).*exp(1).^((1/2).*Di.*l.*v.^(-1))).^2.*(1+exp(1).^((1/2) ...
                        .*Di.*l.*v.^(-1))).^(-2).*Grad.^2.*(exp(1).^((1/2).*Di.*l.*v.^(-1))+(-1) ...
                        .^N.*exp(1).^((1/2).*Di.*(1+N).*l.*v.^(-1)).*((-1)+N)+(-1).^N.*exp(1).^( ...
                        (1/2).*Di.*N.*l.*v.^(-1)).*N).*l.^(-2)+(-2/3).*Di.^(-2).*e.^(-2).*exp( ...
                        1).^((-1/2).*Di.*l.*(2.*e+v.^(-1))).*Grad.^2.*N.*l.^(-2).*v.^(-1).*((-6).* ...
                        ((-1)+exp(1).^(Di.*e.*l)).*((-1).*exp(1).^(Di.*e.*l)+exp(1).^(2.*Di.* ...
                        e.*l)+2.*exp(1).^((1/2).*Di.*l.*v.^(-1))).*v+12.*Di.*e.*exp(1).^((1/2) ...
                        .*Di.*l.*(2.*e+v.^(-1))).*l.*v+Di.^3.*e.^2.*exp(1).^((1/2).*Di.*l.*( ...
                        2.*e+v.^(-1))).*l.^3.*((-3)+8.*e.*v))+(1/96).*Di.^(-2).*e.^(-2).*l.^( ...
                        -2).*v.^(-3).*((-192).*((-1)+exp(1).^((1/4).*Di.*l.*v.^(-1).*(N+(-2).* ...
                        del.*v))).*(e+Grad.^2).*v.^3+(-1).*Di.^3.*l.^3.*(N+(-2).*del.*v).^2.*((-6).* ...
                        e.^2.*v+Grad.^2.*(N+(-2).*del.*v)+e.*(N+4.*del.*v))+(-6).*Di.^2.*l.^2.*v.*(N+ ...
                        (-2).*del.*v).*((-8).*e.^2.*v+Grad.^2.*(N+(-2).*del.*v)+e.*(N+6.*del.*v))+24.* ...
                        exp(1).^((-1).*Di.*del.*l).*v.*((-4).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)).* ...
                        Grad.*v+exp(1).^((1/2).*Di.*del.*l).*Grad.*(4.*v+Di.*l.*(N+(-2).*del.*v))).^2+2.*( ...
                        -1).^(2.*N).*exp(1).^((-1/2).*Di.*del.*l).*Grad.^2.*(48.*exp(1).^((1/4).*Di.* ...
                        N.*l.*v.^(-1)).*v.^2.*((-4).*v+Di.*l.*(N+(-2).*del.*v))+exp(1).^((1/2).* ...
                        Di.*del.*l).*(192.*v.^3+(-6).*Di.^2.*l.^2.*v.*(N+(-2).*del.*v).^2+(-1).* ...
                        Di.^3.*l.^3.*(N+(-2).*del.*v).^3))+(-48).*Di.*exp(1).^((-1/2).*Di.*del.*l).* ...
                        l.*v.^2.*((-4).*e.^2.*(exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).* ...
                        Di.*N.*l.*v.^(-1))).*v+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)).*Grad.^2.*( ...
                        N+(-2).*del.*v)+(-1).*e.*((-4).*exp(1).^((1/2).*Di.*del.*l).*del.*v+exp(1).^( ...
                        (1/4).*Di.*N.*l.*v.^(-1)).*(N+2.*del.*v))));
                    
                    % Summation
                    if any(isnan(W))
                        break
                    else
                        expr = (B./(l.^2)).*W;
                        sumTerm = sumTerm + expr;
                    end
                end
                icsig=exp((-0.5.*gamma.^2./Di.^2).*sumTerm);
                signalModel=...
                    1.*((f.*icsig)+...
                    ((1-f).*exp(-(gamma.^2.*Grad.^2.*(del.^2 .*(DEL-del./3)+...
                    e.^3./30-del.*(e.^2)./6)).*(De./(1+f./2)))));
                %
            case 'DiDe_notort'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                
            case 'DiDe_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
            case 'DiDe_fit_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                fiso=a(6);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*Di);
                    cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
        end
        
        switch ls_mle
            case 'ls'
                % Least squares approach...
                x=sum((sig-signalModel).^2);
            case 'mle'
                besselArg=(signalModel.*sig)./(riceNoiseStd.^2);
                first_part=sum(log(besseli(0,besselArg,1)));
                besselCorrection=sum(abs(real(besselArg)));
                second_part= sum((signalModel.^2)./(2.*(riceNoiseStd.^2)));
                objFun=first_part+besselCorrection-second_part; %function to be maximised
                x= -1.*(objFun); %minimise negative of objFun
        end
    end %End of subFn

%% Signal model formatted for NAG Toolbox e04wd function
    function [mode, x, grad, user]=subFnNag_e04wd(mode, n, a, grad, nstate, user)
        %sumM=0;
        Grad=scanParam(:,1);
        DEL=scanParam(:,2);
        del=scanParam(:,3);
        gamma=scanParam(:,4);
        switch chooseModel
            case 'D'
                R=a(1);
                D=a(2);
                f=a(3);
                S0=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                
            case 'D_fiso'
                R=a(1);
                D=a(2);
                fi=a(3);
                fe=a(4);
                S0=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                
            case 'DiDe'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                
            case 'DiDe_notort'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                
            case 'DiDe_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
            case 'DiDe_fit_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                fiso=a(6);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*Di);
                    cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
        end
        
        % Least squares approach...
        x=sum((sig-signalModel).^2);
    end % End of subFnNag_e04wd

%% Signal model formatted for NAG Toolbox e04us function
    function [mode, f, fjac, user]=subFnNag_e04us(mode, numSubFn, n, ldfj, ...
            needfi, a, fjac, nstate, user)
        %sumM=0;
        Grad=scanParam(:,1);
        DEL=scanParam(:,2);
        del=scanParam(:,3);
        gamma=scanParam(:,4);
        switch chooseModel
            case 'D'
                R=a(1);
                D=a(2);
                f=a(3);
                S0=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                
            case 'D_fiso'
                R=a(1);
                D=a(2);
                fi=a(3);
                fe=a(4);
                S0=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                signalModel=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                
            case 'DiDe'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                
            case 'DiDe_notort'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                
            case 'DiDe_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
            case 'DiDe_fit_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                fiso=a(6);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*Di);
                    cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                signalModel=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
        end
        
        % For e04us, return elements of signalModel, NOT sum of squared
        % differences (with needfi<0, return all elements)
        f=signalModel;
    end % End of subFnNag_e04us

%% As above but for plotting and residuals calculations
    function F=subFnPlot(a,b)
        Grad=b(:,1);
        DEL=b(:,2);
        del=b(:,3);
        gamma=b(:,4);
        switch chooseModel
            case 'D'
                R=a(1);
                D=a(2);
                f=a(3);
                S0=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                F=S0.*((f.*murdayCotts_)+...
                    ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+f./2)))));
                
            case 'D_fiso'
                R=a(1);
                D=a(2);
                fi=a(3);
                fe=a(4);
                S0=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*D);
                    cS=(2+exp(-alpS(m).^2.*D.*(DEL-del)) - 2.*exp(-alpS(m).^2.*D.*del) -2.*exp(-alpS(m).^2.*D.*DEL) + exp(-alpS(m).^2.*D.*(DEL+del)) )./((alpS(m).^2.*D).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*D)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*D.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*D.*del_vec) -2.*exp(-alpS_vsq.*D.*DEL_vec) + exp(-alpS_vsq.*D.*(DEL_vec+del_vec)) )./((alpS_vsq.*D).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                
                F=S0.*( (fi.*murdayCotts_) + ...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(D))));
                
            case {'DiDe','DiDe_unnormalised'}
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                switch chooseModel
                    case 'DiDe'
                        % nothing to do
                    case 'DiDe_unnormalised'
                        s0=a(5);
                end
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                switch chooseModel
                    case 'DiDe'
                        F=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                    case 'DiDe_unnormalised'
                        F=s0.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                end
                
            case 'DiDe_with_rise_time'
                %
                e=scanParam(:,5);
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                sumTerm = 0;
                N = 1; % Number of half oscillations (i.e. N = 1 for PGSE)
                v = N./(2.*del);
                for n=1:numel(rootsS)
                    u = rootsS(n);
                    l = (u/R)^2;
                    B = (2*(R/u)^2)/(u^2 - 2);
                    % Waveform
                    % Gamma2 in GPD_trapezoid_arbitrary_freq.m - using ToMatlab package
                    W=(-2).*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^(-1).*(N+2.*(del+DEL).* ...
                        v)).*(exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)) ...
                        ).^4.*Grad.^2.*l.^(-2)+2.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^( ...
                        -1).*(N+2.*(2.*e+DEL).*v)).*((-1)+exp(1).^(Di.*e.*l)).^2.*(exp(1).^(Di.* ...
                        e.*l)+(-1).*exp(1).^((1/2).*Di.*l.*v.^(-1))).^2.*(1+exp(1).^((1/2).* ...
                        Di.*l.*v.^(-1))).^(-2).*((-1).^(1+N)+exp(1).^((1/2).*Di.*N.*l.*v.^(-1))) ...
                        .*((-1)+(-1).^N.*exp(1).^((1/2).*Di.*N.*l.*v.^(-1))).*Grad.^2.*l.^(-2)+2.*( ...
                        -1).^N.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.*l.*v.^(-1).*(N+2.*(e+ ...
                        del+DEL).*v)).*((-1)+exp(1).^(Di.*e.*l)).*(exp(1).^(Di.*e.*l)+(-1).*exp(1) ...
                        .^((1/2).*Di.*l.*v.^(-1))).*(1+exp(1).^((1/2).*Di.*l.*v.^(-1))).^(-1).*( ...
                        exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1))).^2.* ...
                        (1+(-1).^(1+N).*exp(1).^(Di.*del.*l)+(-2).*exp(1).^(Di.*DEL.*l)+exp(1).^( ...
                        Di.*del.*l+(1/2).*Di.*N.*l.*v.^(-1))+2.*(-1).^N.*exp(1).^(Di.*DEL.*l+(1/2).* ...
                        Di.*N.*l.*v.^(-1))+(-1).^(1+N).*exp(1).^((1/2).*Di.*N.*l.*v.^(-1))).* ...
                        Grad.^2.*l.^(-2)+(-4).*(-1).^N.*Di.^(-2).*e.^(-2).*exp(1).^((-1/2).*Di.* ...
                        l.*v.^(-1).*(1+N+4.*e.*v)).*((-1)+exp(1).^(Di.*e.*l)).^2.*(exp(1).^( ...
                        Di.*e.*l)+(-1).*exp(1).^((1/2).*Di.*l.*v.^(-1))).^2.*(1+exp(1).^((1/2) ...
                        .*Di.*l.*v.^(-1))).^(-2).*Grad.^2.*(exp(1).^((1/2).*Di.*l.*v.^(-1))+(-1) ...
                        .^N.*exp(1).^((1/2).*Di.*(1+N).*l.*v.^(-1)).*((-1)+N)+(-1).^N.*exp(1).^( ...
                        (1/2).*Di.*N.*l.*v.^(-1)).*N).*l.^(-2)+(-2/3).*Di.^(-2).*e.^(-2).*exp( ...
                        1).^((-1/2).*Di.*l.*(2.*e+v.^(-1))).*Grad.^2.*N.*l.^(-2).*v.^(-1).*((-6).* ...
                        ((-1)+exp(1).^(Di.*e.*l)).*((-1).*exp(1).^(Di.*e.*l)+exp(1).^(2.*Di.* ...
                        e.*l)+2.*exp(1).^((1/2).*Di.*l.*v.^(-1))).*v+12.*Di.*e.*exp(1).^((1/2) ...
                        .*Di.*l.*(2.*e+v.^(-1))).*l.*v+Di.^3.*e.^2.*exp(1).^((1/2).*Di.*l.*( ...
                        2.*e+v.^(-1))).*l.^3.*((-3)+8.*e.*v))+(1/96).*Di.^(-2).*e.^(-2).*l.^( ...
                        -2).*v.^(-3).*((-192).*((-1)+exp(1).^((1/4).*Di.*l.*v.^(-1).*(N+(-2).* ...
                        del.*v))).*(e+Grad.^2).*v.^3+(-1).*Di.^3.*l.^3.*(N+(-2).*del.*v).^2.*((-6).* ...
                        e.^2.*v+Grad.^2.*(N+(-2).*del.*v)+e.*(N+4.*del.*v))+(-6).*Di.^2.*l.^2.*v.*(N+ ...
                        (-2).*del.*v).*((-8).*e.^2.*v+Grad.^2.*(N+(-2).*del.*v)+e.*(N+6.*del.*v))+24.* ...
                        exp(1).^((-1).*Di.*del.*l).*v.*((-4).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)).* ...
                        Grad.*v+exp(1).^((1/2).*Di.*del.*l).*Grad.*(4.*v+Di.*l.*(N+(-2).*del.*v))).^2+2.*( ...
                        -1).^(2.*N).*exp(1).^((-1/2).*Di.*del.*l).*Grad.^2.*(48.*exp(1).^((1/4).*Di.* ...
                        N.*l.*v.^(-1)).*v.^2.*((-4).*v+Di.*l.*(N+(-2).*del.*v))+exp(1).^((1/2).* ...
                        Di.*del.*l).*(192.*v.^3+(-6).*Di.^2.*l.^2.*v.*(N+(-2).*del.*v).^2+(-1).* ...
                        Di.^3.*l.^3.*(N+(-2).*del.*v).^3))+(-48).*Di.*exp(1).^((-1/2).*Di.*del.*l).* ...
                        l.*v.^2.*((-4).*e.^2.*(exp(1).^((1/2).*Di.*del.*l)+(-1).*exp(1).^((1/4).* ...
                        Di.*N.*l.*v.^(-1))).*v+(-1).*exp(1).^((1/4).*Di.*N.*l.*v.^(-1)).*Grad.^2.*( ...
                        N+(-2).*del.*v)+(-1).*e.*((-4).*exp(1).^((1/2).*Di.*del.*l).*del.*v+exp(1).^( ...
                        (1/4).*Di.*N.*l.*v.^(-1)).*(N+2.*del.*v))));
                    
                    % Summation
                    if any(isnan(W))
                        break
                    else
                        expr = (B./(l.^2)).*W;
                        sumTerm = sumTerm + expr;
                    end
                end
                icsig=exp((-0.5.*gamma.^2./Di.^2).*sumTerm);
                F=...
                    1.*((f.*icsig)+...
                    ((1-f).*exp(-(gamma.^2.*Grad.^2.*(del.^2 .*(DEL-del./3)+...
                    e.^3./30-del.*(e.^2)./6)).*(De./(1+f./2)))));
                
            case 'DiDe_notort'
                R=a(1);
                Di=a(2);
                De=a(3);
                f=a(4);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts_=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=1.*((f.*murdayCotts_)+...
                ((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                F=1.*((f.*murdayCotts)+((1-f).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                
            case 'DiDe_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                bS=(2.*del)./(alpS(m).^2.*Di);
                cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                exprM=aS.*(bS-cS);
                sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=1.*((fi.*murdayCotts)+...
                (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                F=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((1-fi-fe).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
            case 'DiDe_fit_fiso'
                R=a(1);
                Di=a(2);
                De=a(3);
                fi=a(4);
                fe=a(5);
                fiso=a(6);
                alpS=rootsS./R;
                %{
                sumM=0;
                for m=1:numel(rootsS)
                    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
                    bS=(2.*del)./(alpS(m).^2.*Di);
                    cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
                    exprM=aS.*(bS-cS);
                    sumM=sumM+exprM;
                end
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumM);
                F=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
                %}
                nroot=numel(rootsS);
                ndel=numel(del);
                
                % Vectorize assignments
                aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
                bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);
                
                % Duplicate vectors into matrices to enable vectorised assignment
                alpS_vsq=repmat(alpS,1,ndel)'.^2;
                del_vec=repmat(del,1,nroot);
                DEL_vec=repmat(DEL,1,nroot);
                cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);
                
                aSvec_vec=repmat(aSvec,1,ndel)';
                exprM_vec=aSvec_vec.*(bSvec-cSvec);
                
                sumMvec = sum(exprM_vec,2);
                murdayCotts=exp(-2.*(gamma.^2).*(Grad.^2).*sumMvec);
                F=1.*((fi.*murdayCotts)+...
                    (fe.*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+fi./2))))+...
                    ((fiso).*exp(-(((Grad.*del.*gamma).^2).*(DEL-del./3)).*(De))));
        end
    end % End of subFnPlot

end % End of main function