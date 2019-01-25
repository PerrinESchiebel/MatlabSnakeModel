%%edited from MultiPeg_Flexible_v3--STARTED 082017 BY Perrin
%%Switch to constant Do, Ro, sweep different kplds
%%%08/24/17  fixed wavephase. Centered initial conditions box around zero.
%%%Removed no return condition as some trajectories were re-entering pegs
%%%09/07/17 added to find all points of contact with post
%%%09/12/17 changed to find Po every time it touches a peg
%%%09/14/17 changed kappa_max to match winding number of 6.2
%%%01/15/18 Realized that the wavenumber vs length vs peg radius matters
%%%since I'm using real units, esp. in terms of comparing to Zeb Rocklin's
%%%models. Changed length to 0.2 to match wavenumber of 1 but I haven't
%%%re-run anything. I can maybe change peg radius. I should also look at
%%%whether peg radius really matters, or if it's more distance between
%%%them/initial condition that sets scattering
% DoLIST = [0.07,0.08];
addpath F:\Dropbox\Research\scattering\snake_mech
directory = pwd;
hh = strfind(directory,'Dropbox');
directory = directory(1:hh-1);

DoLIST = [0.022,0.023,0.025, 0.05];
DoLIST = 0.023;
% DoLIST = [0.025,0.027,0.03,0.035,0.04,0.05];
% DoLIST = 0.023;
% % DoLIST = 0.05;
% % Ro = 0.015;
% % RoLIST = 0.0032;
% % RoLIST = 0.0082;
% % RoLIST = [0.0042];
% % RoLIST = 0.01;
% % RoLIST = 0.05;
% % DoLIST = 0.3;
RoLIST = [0.01];
% % RoLIST = 0.005;
% % RoLIST = 0.00032;
% % RoLIST = 0.0038;
% DoLIST = 0.023;
% DoLIST = 0.023/2;
% RoLIST = 0.01/2;
% % DoLIST = 0.0154-0.00032*2;
% % RoLIST = (1/1000).*(23-15.4)/2;
% % DoLIST = 0.023;
% RoLIST = 0.0254;
% kpldsLIST = 9.3;
% kpldsLIST = [3,7];
kpldsLIST = 5;
kpldMAXLIST = [5, 10, 40, 150, 500];
% DoLIST = 0.075;
% RoLIST = 0.026;
%                             plotyesno = true;
str = input('Do you want to save these runs?','s');
switch str
    case 'y'
        saveyesno = true;
    case 'n'
        saveyesno = false;
end
%%
% DoLIST = [0.0525,0.057,0.061,0.069,0.076] + Ro;
for mm = 1:length(DoLIST)
    for ww = 1:length(kpldsLIST)
        for rr = 1:length(RoLIST)
            %%
            tic
%                                     plotyesno = true;
            plotyesno = false;
            repeat = 7;
            % repeat = 3;
            npegs = 9;
            %         touchradius = 0.003;
            touchradius = 0.003;
            %%robot  shape
            % L = 0.8;
            % % vc = -0.12;  %%for frequency of 0.15 Hz, vc = f*L/xi = 0.15*0.8/1
            % vc = -1;
            % thm = 9.2/(2*pi);
            % xi = 1;
            %%snake shape
%                         L = 0.4; %%realized that L and xi set the size of the wave, I
            %             can use xi = 1 but need to adjust length by 2
            L = 0.2;
            vc = -1;
            thm = kpldsLIST(ww)/(2*pi);  %%kappalambda = 5
            xi = 1;
            
            %%%%
            lambdas = L/xi;
            periodinframes = 500;
            period = abs(L/(xi*vc));
            dt = period/periodinframes;
            ds = abs(vc*dt);
            s = 0:ds:L;
            numsegs = length(s);
            numtimes = repeat*periodinframes;
            t = linspace(dt,repeat*period,numtimes);%%%%this needs to start at dt NOT 0
            %%%%CALCULATE THE KINEMATICS FOR THE HOMOG SUBSTRATE CASE
            [ss, tt] = meshgrid(s, t);
            th = thm.*sin(2 * pi / lambdas * (ss - vc * tt));
            tmpz = cumtrapz(s, cos(th'))';
            tmpx = cumtrapz(s, sin(th'))';
            zloc_mean = mean(tmpz, 2);
            xloc_mean = mean(tmpx, 2);
            z = tmpz - repmat(zloc_mean, [1, numsegs]);
            x = tmpx - repmat(xloc_mean, [1, numsegs]);
            
            th0 = thm.*sin(2 * pi / lambdas * (zeros(numtimes, numsegs) - vc * tt));
            vzLoc = -vc .* cos(th) + vc .* cos(th0);
            vxLoc = -vc .* sin(th) + vc .* sin(th0);
            vz =  vzLoc - repmat(mean(vzLoc, 2), [1 numsegs]);
            vx =  vxLoc - repmat(mean(vxLoc, 2), [1 numsegs]);
            dzetadt = -thm*2*pi/lambdas*2*pi/period.*sin(2*pi/lambdas*(ss-vc*tt)).*ds;
            %            shape = [thm,2,L,1/period];
            %%
            %                         [F,Fnormal,~] = RFT(thm,2,L,1/period,round(numsegs*2/xi),periodinframes); %%% need to adjust number of segments so that first numsegs amount corresponds to xi = 1
            % %             F = repmat(F,repeat,1);
            %             Fnormal = repmat(Fnormal,repeat,1);
            %             Fnormal = Fnormal(:,1:numsegs,:);
            % %             F = F(:,1:numsegs,:);
            %%%%THIS KAPPA IS THE CONTROL TARGET FOR THE SNAKE AT ALL TIME
            %%%NOTE KAPPA*DS = ZETA
            kappadesired = thm*2*pi/lambdas.*cos(2 * pi / lambdas * (ss - vc * tt));
            %%%%CoM VELOCITY THAT MINIMIZES SLIP
            vz = vz(1,:);vx = vx(1,:);
            vcm = 1-min(sqrt(vz(1,:).^2+vx(1,:).^2));  %%com velocity for no-slip
            for jj=1:numtimes
                z(jj,:) = z(jj,:) + vcm*dt*jj;
            end
            %%
            %%%%%%%%%%%%%%%%%%%SCATTERING INTERACTION%%%%%%%%%%%%%%%%%
            %%%CREATE INITIAL CONDITIONS BOX%%%%%%
            %%%%INITIAL CONDITIONS BOX ACHIEVED BY MOVING PEG
            Do = DoLIST(mm);
            %             disp(DoLIST(mm));
            Ro = RoLIST(rr);
            radii = repmat(Ro,npegs,1);
            %             nzs = 100;
%             nzs = 100;
            nzs = 100;
            nsegs1wave = round(numsegs/xi);
            deltaz = z(1,nsegs1wave) - z(1,1);
            deltax = Do;
            initz = linspace(deltaz/2+Ro+0.01,0.01+deltaz/2+Ro+deltaz,nzs);
            %         nxs = ceil(deltax/mean(diff(initz)));
            dx = mean(diff(initz));
            %         initx = (-deltax/2):dx:(deltax/2);
            %         initx = 0:dx:deltax/2;
            initx = [-fliplr(dx:dx:deltax/2),0:dx:deltax/2];
            nxs = length(initx);
%             initx = 0;
%             nxs = 1;
            %         initx = linspace(-deltax/2,deltax/2,nxs);
            %             zeta_max_pervertebrae = 0.105*[-1;1]; %%rads per "vertebrae". Assume there are 100 vertebrae involved in one wave.
            %             kappa_max = 100/(xi*numsegs)*zeta_max_pervertebrae/ds;  %%curvature per segment given xi*numsegs/100 segments per vertebrae
%             zeta_max_pervertebrae = 0.105*[-1;1]; %%rads per "vertebrae". Assume there are 100 vertebrae involved in one wave.
%             kappa_max = 93/(xi*numsegs)*zeta_max_pervertebrae/ds;  %%curvature per segment given xi*numsegs/100 segments per vertebrae
%            kappa_max = 80*[-1;1];
            radsperm = 101;radsperbody = radsperm*L;
            radspersegment = radsperbody/numsegs;
            kappa_max = (radspersegment/ds)*[-1;1];   %%%using kappa = dtheta/ds   and radspersegment = dtheta
            kappa_max = 50;
            x0 = [0.8 0];
            count = 0;
            numruns = nxs*nzs;
            savehead = struct('pegR',[],'pegpos',[],'x',[],'y',[]);
            pegphaseALL = cell(numruns,1);
            pegorderdata = nan(numruns,npegs);
            alltimetouch = false(numruns,numtimes);
            wavephase = nan(numruns,1);pegphase = nan(numruns,1);countercwORcw = false(numruns,2);
            AllBodyThetaPeg = cell(numruns,1);
            AllsegmentTouches = cell(numruns,1);
            thetadifferenceALL = nan(numruns,1);
            options = optimset('Display','Off');
%             for ll = 1:nxs
                            for ll = 1
%                 for kk = 1:nzs
                                    for kk=10
                    disp(count);
                    count = count+1;
                    spacings = (-(npegs-1)/2:(npegs-1)/2)';
                    pegpositions= [repmat(initz(kk),npegs,1),repmat(initx(ll),npegs,1)+spacings*Do];
                    znext = z(1,end);                 %%%SET BEGINNING HEAD POSITION
                    xnext = x(1,end);
                    thpred = th(1,:);
                    zact = nan(size(z));            %%%SAVE THE CALCULATED TRAJECTORY IN XACT,YACT
                    xact = nan(size(z));
                    kapparealized = nan(size(th));
                    kappapred = kappadesired(1,:);
                    thetaactual = nan(size(th));
                    vzsave = zeros(1,numsegs);
                    vxsave = zeros(1,numsegs);
                    touchedApegEVER = false;
                    firstTouch = true;
                    touchingnow = false;
                    pegtouch = false(1,npegs);
                    pegtouchorder = zeros(1,npegs);
                    allthetapegs = zeros(numtimes,3);pind = 1;
                    timefirsttouch = 0;pegphasefirsttouch = 0;
                    CCWorCWperpeg = false(npegs,2);
                    ABTP = cell(numtimes,1);
                    segmentsthattouched = false(numtimes,npegs,numsegs);
                    for jj = 1:numtimes-1
                        pegIn = sqrt((pegpositions(:,1) - znext).^2 + (pegpositions(:,2) - xnext).^2) < Ro;
                        if any(pegIn)
                            touchedApegEVER = true;
                            touchedApegNOW = true;
                        else
                            touchedApegNOW = false;
                        end
                        if touchedApegEVER == false && firstTouch == true
                            znext = z(jj+1,end);xnext = x(jj+1,end);
                            xact(jj,:) = x(jj,:);zact(jj,:) = z(jj,:);thact = th(jj,:);
                            thpred = th(jj+1,:);
                            kappapred = kappadesired(jj+1,:);
                            vzsave(jj) = vcm;
                            thetaactual(jj,:) = thact;
                            if plotyesno
                                %                             figure(2);plot(zact(jj,:),xact(jj,:),'LineWidth',3);hold on;viscircles(pegpositions,radii);
                                %                             hold off;drawnow;axis equal tight;
                            end
                        else
                            if touchedApegNOW == true && firstTouch == true
                                touchtime1 = jj;
                                getv = true;
                                touchingnow = true;
                                %%FIND WHERE SNAKE HEAD CONTACTS PEG
                                xt = pegpositions(pegIn,2) - xnext;
                                zt =  pegpositions(pegIn,1) - znext;
                                pegtouch(pegIn) = 1;
                                peglocationCurrent = pegpositions(pegIn,:);
                                pegtouchorder = pegtouchorder + pegtouch;
                                thetapeg(1) = atan2(-zt,xt);
                                thetapeg(2) = atan2(zt,-xt);
                                thetadifference = thetapeg-thpred(end);
                                %%THIS IS THE DIRECTION THE SNAKE WILL GO AROUND THE PEG
                                CCWorCW = abs(thetadifference) == min(abs(thetadifference));
                                thetapeg = thetapeg(CCWorCW);
                                thetadifferenceZR = thetadifference(CCWorCW);
                                %%%%%SAVE INFO FOR PHASES
                                CCWorCWperpeg(pegIn,:) = CCWorCW;
                                if timefirsttouch == 0 && pegphasefirsttouch == 0
                                    timefirsttouch = jj;
                                    pegphasefirsttouch = thetapeg;
                                    countercwORcw(count,:) = CCWorCW;
                                end
                                %%FIND Po
                                %                                 if CCWorCW(1)
                                %                                     Po = find(thpred>0,1,'last');
                                %                                 elseif CCWorCW(2)
                                %                                     Po = find(thpred<0,1,'last');
                                %                                 end
                                firstTouch = false;
                            elseif touchedApegNOW == true && firstTouch == false
                                xt = peglocationCurrent(2) - xnext;
                                zt =  peglocationCurrent(1) - znext;
                                thetapeg(1) = atan2(-zt,xt);
                                thetapeg(2) = atan2(zt,-xt);
                                thetapeg = thetapeg(CCWorCW);
                                thetadifference = CCWorCW*[-1;1]*(thetapeg-thpred(end));
                                if thetadifference <= 0
                                    firstTouch = true;
                                    touchingnow = false;
                                    %                                 noReturns(pegtouch) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %                             display('left peg');
                                end
                            end
                            if touchingnow == true && touchedApegNOW == true
                                alltimetouch(count,jj) = true;
                                %                                 allsegtouch(count,jj) = true;
                                kappaedit = kappapred;
                                if CCWorCW(1)
                                    Po = find(thpred>0,1,'last');
                                elseif CCWorCW(2)
                                    Po = find(thpred<0,1,'last');
                                end
                                kappaedit(Po) = CCWorCW*kappa_max;
                                thact = cumtrapz(s,kappaedit) + th0(jj,:);
                                tmpz = cumtrapz(s, cos(thact'))';
                                tmpx = cumtrapz(s, sin(thact'))';
                                zloc_mean = mean(tmpz(:,1:nsegs1wave), 2);
                                xloc_mean = mean(tmpx(:,1:nsegs1wave), 2);
                                ztest = tmpz - repmat(zloc_mean, [1,numsegs]);
                                xtest = tmpx - repmat(xloc_mean, [1,numsegs]);
                                vzLoc = -vc .* cos(thact) + vc .* cos(th0(jj,:));
                                vxLoc = -vc .* sin(thact) + vc .* sin(th0(jj,:));
                                vz =  vzLoc - repmat(mean(vzLoc, 2), [1 numsegs]);
                                vx =  vxLoc - repmat(mean(vxLoc, 2), [1 numsegs]);
                                %%%%FIND THE CoM V THAT MINIMIZES SLIP
                                funv = @(V) 1 - sum(dotProducts(V,vz,vx,cos(thact),sin(thact)))./numsegs;
                                [v0] = fminsearch(funv,x0);
                                x0 = v0;
                                vzsave(jj) = v0(1);
                                vxsave(jj) = v0(2);
                                ztest = ztest + sum(vzsave.*dt);
                                xtest = xtest + sum(vxsave.*dt);
                                zhead = ztest(end);xhead = xtest(end);
                                %%SEE IF THE SNAKE HEAD IS STILL IN THE PEG
                                xt = peglocationCurrent(2) - xhead;
                                zt =  peglocationCurrent(1) - zhead;
                                thetapeg(1) = atan2(-zt,xt);
                                thetapeg(2) = atan2(zt,-xt);
                                thetapeg = thetapeg(CCWorCW);
                                thetadifference = CCWorCW*[-1;1]*(thetapeg-thact(end)) >0;
                                pegdis = sqrt((peglocationCurrent(1) - zhead).^2 + (peglocationCurrent(2) - xhead).^2) < Ro+Ro*0.1;
                                nsegs = 0;
                                while thetadifference && pegdis
                                    nsegs = nsegs+1;
                                    clear xtest ztest
                                    kappaedit(Po-nsegs) = CCWorCW*kappa_max;
                                    thact = cumtrapz(s,kappaedit) + th0(jj,:);
                                    tmpz = cumtrapz(s, cos(thact'))';
                                    tmpx = cumtrapz(s, sin(thact'))';
                                    zloc_mean = mean(tmpz(:,1:nsegs1wave), 2);
                                    xloc_mean = mean(tmpx(:,1:nsegs1wave), 2);
                                    ztest = tmpz - repmat(zloc_mean, [1,numsegs]);
                                    xtest = tmpx - repmat(xloc_mean, [1,numsegs]);
                                    vzLoc = -vc .* cos(thact) + vc .* cos(th0(jj,:));
                                    vxLoc = -vc .* sin(thact) + vc .* sin(th0(jj,:));
                                    vz =  vzLoc - repmat(mean(vzLoc, 2), [1 numsegs]);
                                    vx =  vxLoc - repmat(mean(vxLoc, 2), [1 numsegs]);
                                    %%%%FIND THE CoM V THAT MINIMIZES SLIP
                                    funv = @(V) 1 - sum(dotProducts(V,vz,vx,cos(thact),sin(thact)))./numsegs;
                                    [v0] = fminsearch(funv,x0);
                                    x0 = v0;
                                    vzsave(jj) = v0(1);
                                    vxsave(jj) = v0(2);
                                    ztest = ztest + sum(vzsave.*dt);
                                    xtest = xtest + sum(vxsave.*dt);
                                    zhead = ztest(end);xhead = xtest(end);
                                    xt = peglocationCurrent(2) - xhead;
                                    zt =  peglocationCurrent(1) - zhead;
                                    thetapeg(1) = atan2(-zt,xt);
                                    thetapeg(2) = atan2(zt,-xt);
                                    thetapeg = thetapeg(CCWorCW);
                                    thetadifference = CCWorCW*[-1;1]*(thetapeg-thact(end)) >0;
                                    pegdis = sqrt((peglocationCurrent(1) - zhead).^2 + (peglocationCurrent(2) - xhead).^2) < Ro;
                                end
                                allthetapegs(pind,1) = find(pegIn == true);
                                allthetapegs(pind,2) = jj;
                                allthetapegs(pind,3) = thetapeg;
                                pind = pind+1;
                                zact(jj,:) = ztest;
                                xact(jj,:) = xtest;
                                kapparealized(jj,:) = kappaedit;
                                kappapred(1:end-1) = kappaedit(2:end);
                                kappapred(end) = kappadesired(jj+1,end);
                                thetaactual(jj,:) = thact;
                                %%
                                %%%%SET THE PREDICTION FOR THE HEAD TO THE NEXT THETA_DESIRED
                                %%%% CALCULATE THE NEW SHAPE GIVEN PREDICTED THETA SO WE
                                %%%% CAN CHECK IF IT WILL BE INSIDE THE PEG AGAIN
                                thpred = nan(1,length(thact));
                                thpred(1:end-1) = thact(2:end);
                                thpred(end) = kappadesired(jj+1,end)*ds + thact(end);
                                tmpz = cumtrapz(s, cos(thact'))';
                                tmpx = cumtrapz(s, sin(thact'))';
                                zloc_mean = mean(tmpz(:,1:nsegs1wave), 2);
                                xloc_mean = mean(tmpx(:,1:nsegs1wave), 2);
                                zpred = tmpz - repmat(zloc_mean, [1, numsegs]);
                                xpred = tmpx - repmat(xloc_mean, [1, numsegs]);
                                zpred = zpred + sum(vzsave.*dt);
                                xpred = xpred + sum(vxsave.*dt);
                                znext = zpred(end);xnext = xpred(end);
                                if plotyesno
                                    figure(2);plot(zact(jj,:),xact(jj,:),'k','LineWidth',3);hold on;viscircles(pegpositions,radii);
                                    %                                                 quiver(zact(jj,end),xact(jj,end),cos(thetapeg),sin(thetapeg),0.1,'Color','r','LineWidth',2);
                                    %                     quiver(zact(jj,end),xact(jj,end),cos(thetapeg(2)),sin(thetapeg(2)),0.1,'Color','c','LineWidth',2);
                                    quiver(zact(jj,end),xact(jj,end),cos(thact(end)),sin(thact(end)),0.1,'Color','k','LineWidth',2);
                                    %                                                 quiver(zact(jj,end),xact(jj,end),cos(thpred(end)),sin(thpred(end)),0.1,'Color','g','LineWidth',2);
                                    viscircles(pegpositions(pegtouch,:),radii(pegtouch),'EdgeColor',[0 0.5 0.5]);
                                    %                                                 plot(zact(jj,:),xact(jj,:),'LineWidth',3,'Color',[0.2 0.2 0.2]);
                                    axis equal tight;drawnow;hold off;pause(0.5);
                                end
                            elseif touchingnow == false  || touchedApegNOW == false    %%%IF THE SNAKE ALREADY CONTACTED THE PEG BUT ISN'T ANYMORE
                                kapparealized(jj,:) = kappapred;
                                thact = thpred; %%%KEEP PASSING THE CURRENT SHAPE DOWN
                                tmpz = cumtrapz(s, cos(thact'))';
                                tmpx = cumtrapz(s, sin(thact'))';
                                zloc_mean = mean(tmpz(:,1:nsegs1wave), 2);
                                xloc_mean = mean(tmpx(:,1:nsegs1wave), 2);
                                zact(jj,:) = tmpz - repmat(zloc_mean, [1,numsegs]);
                                xact(jj,:) = tmpx - repmat(xloc_mean, [1,numsegs]);
                                vzLoc = -vc .* cos(thact) + vc .* cos(th0(jj,:));
                                vxLoc = -vc .* sin(thact) + vc .* sin(th0(jj,:));
                                vz =  vzLoc - repmat(mean(vzLoc, 2), [1 numsegs]);
                                vx =  vxLoc - repmat(mean(vxLoc, 2), [1 numsegs]);
                                if getv
                                    %                             display(getv)
                                    funv = @(V) 1 - sum(dotProducts(V,vz,vx,cos(thact),sin(thact)))./numsegs;
                                    [v0,~,exitflag] = fminsearch(funv,x0,options);
                                    if  jj - touchtime1 >numsegs+100
                                        getv = false;
                                        display(getv);
                                        v0(1) = mean(vzsave(jj-100:jj-1));
                                        v0(2) = mean(vxsave(jj-100:jj-1));
                                    end
                                    x0 = v0;
                                end
                                vzsave(jj) = v0(1);
                                vxsave(jj) = v0(2);
                                zact(jj,:) = zact(jj,:) + sum(vzsave.*dt);
                                xact(jj,:) = xact(jj,:) + sum(vxsave.*dt);
                                znext = zact(jj,end);xnext = xact(jj,end);
                                thetaactual(jj,:) = thact;
                                kappapred(1:end-1) = kappapred(2:end);
                                kappapred(end) = kappadesired(jj+1,end); %%HAVE THE HEAD GO ALONG ON NOMINAL SHAPE
                                thpred(1:end-1) = thpred(2:end); %%%KEEP PASSING THE CURRENT SHAPE DOWN
                                thpred(end) = kappadesired(jj+1,end)*ds + thpred(end-1);
                                if plotyesno
                                    figure(2);plot(zact(jj,:),xact(jj,:),'c','LineWidth',3);hold on;viscircles(pegpositions,radii);
                                    viscircles(pegpositions(pegtouch,:),radii(pegtouch),'EdgeColor',[0.1 0.5 0.5]);axis equal tight;drawnow;hold off;pause(0.5);
                                end
                            end
                        end
                        clear anypegIn
                        anypegIn = sqrt((pegpositions(:,1) - zact(jj,:)).^2 + (pegpositions(:,2) - xact(jj,:)).^2) < (Ro+touchradius);
                        if touchedApegEVER
                            if any(any(anypegIn'))
                                postsTouched = find((any(anypegIn')).*pegtouch);
                                ABTPnow = cell(npegs,1);
                                
                                for pt = postsTouched
                                    peglocationCurrent = pegpositions(pt,:);
                                    xt = peglocationCurrent(2) - xact(jj,anypegIn(pt,:));
                                    zt =  peglocationCurrent(1) - zact(jj,anypegIn(pt,:));
                                    thetapegccw = atan2(-zt,xt);
                                    thetapegcw = atan2(zt,-xt);
                                    thetapegtest = [thetapegccw;thetapegcw];
                                    oneortwo = find(CCWorCWperpeg(pt,:));
                                    thetapegtest = thetapegtest(oneortwo,:);
                                    ABTPnow{pt} = thetapegtest';
                                    segmentsthattouched(jj,pt,:) = anypegIn(pt,:);
                                    if jj>50
                                        
                                    end
                                end
                                ABTP{jj} = cell2mat(ABTPnow);
                                
                            end
                        end
                    end
                    savehead(count).pegpos = [initz(kk),initx(ll)];
                    savehead(count).pegR = Ro;
                    savehead(count).z = zact(:,1);
                    savehead(count).x = xact(:,1);
                    [~,ind] = sort((pegtouchorder),'descend');
                    pegorderdata(count,:) = ind.*pegtouch(ind);
                    %                 wavephase(count) = mod(-timefirsttouch/periodinframes * vc * 2 * pi / lambdas,2*pi); %%
                    wavephase(count) = -timefirsttouch/periodinframes; %%
                    pegphase(count) = pegphasefirsttouch;
                    pegphaseALL{count} = allthetapegs(allthetapegs(:,1)>0,:);
                    AllBodyThetaPeg{count} = cell2mat(ABTP);
                    AllsegmentTouches{count} = segmentsthattouched;
                    thetadifferenceALL(count) = thetadifferenceZR;
                end
            end
            zz = nan(length(savehead),numtimes);
            for jj=1:length(savehead)
                zz(jj,:) = savehead(jj).z - savehead(jj).pegpos(1);
            end
            xx = nan(length(savehead),numtimes);
            for jj=1:length(savehead)
                xx(jj,:) = savehead(jj).x-savehead(jj).pegpos(2);
            end
            FindScatteringAngle
            if saveyesno
                save([directory,'Dropbox\Research\scattering\snake_mech\results\DoSweep_incrkappamax\','kappam_',num2str(round(kappa_max(2))),'_L0p2',num2str(npegs),'_Do_',num2str(DoLIST(mm)*1000),'_Ro_',num2str(Ro*1000),'_kpl_',num2str(kpldsLIST(ww)),'.mat']);
                clearvars -except DoLIST kpldsLIST directory RoLIST mm ww rr saveyesno
            end
            toc
        end
    end
end