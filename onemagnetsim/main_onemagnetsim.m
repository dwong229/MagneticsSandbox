%% Model 1 magnet and 1 electromagnetic coil in 1D

close all
clear all

simMode = 3;
graphicalsimualtion = true;
% System definition
% -----------------------
%% Sliding friction model
% mu_static = 0.005;
% mu_kinetic = 0.001;
% g = 9.81; %acceleration due to gravity
% N = 0.1*g;
% F_drag = mu_kinetic*N*sign(v)*(-1);

global rho magnet coil 

%% Viscous drag friction model
rho = 1000; % density of water 1000kg/m3
mass = 0.1;
R = 0.01;
Cd = 1.05*10;
u1 = 1;
%package variables
magnet.mass = mass; %mass of the magnet, kg
magnet.R = 0.01;
magnet.A = magnet.R^2*pi;
%Cd = 0.47; % drag coeff for sphere
magnet.Cd = Cd; % drag coeff for square
magnet.ui = u1;% magnetic dipole moment of permanet magnet (Am^2)

% electromagnetic coil parameters
u0 = 4*pi*10^-7;
R1 = 0.1;
I1 = 50;

coil.u0 = 4*pi*10^-7; % permiability of free space
coil.R = R1; %coil radius m
coil.I = I1; % apply a constant current
coil.x = 0;

switch simMode
    case 1
        disp('1 coil 1 magnet - constant current')        
        % initial condition
        % ----------------------
        x0_mag = 0.5;
        v0 = 0;
        
        X = [x0_mag v0];
        
        % Timing
        % ------
        tStart = 0;  % simulation start time (sec)
        tFinal = 100; % simulation start time (sec)
        
        
        graphicalsimualtion = true;
        %% Physics simulator
        % Set up the magnets starting config;
        Xinitial = X(end,:);
        
        tspan = [tStart tFinal];
        
        % solve for posn using ode45
        %[tvector, Xmatrix, timeAtEvent, stateAtEvent] = ode45(@compute_magnet_derivatives, tspan, Xinitial, option);
        [tvector, Xmatrix] = ode45(@compute_magnet_derivatives, tspan, Xinitial);
        sol = ode45(@compute_magnet_derivatives, tspan, Xinitial);
        
        [xxdot,xdotxddot] = deval(sol,tvector);
        
        force = xdotxddot(2,:)/magnet.mass;
        
        figure
        subplot(3,1,1)
        plot(tvector,Xmatrix(:,1));
        title('Position')
        subplot(3,1,2)
        plot(tvector,Xmatrix(:,2));
        title('Velocity')
        
        subplot(3,1,3)
        plot(tvector,xdotxddot(2,:));
        title('Acceleration')
        xlabel('time (sec)')
        
        %% Visualization
        
        if graphicalsimualtion
            
            h1 = figure;
            hcoil = rectangle('Position',[coil.x-coil.R/8,-coil.R,coil.R/4,coil.R*2],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[coil.I/coil.I 0 0]);
            hold on
            axis equal
            
            hmagnet = drawmagnet(Xmatrix(1,1),0,magnet.R*4,magnet.R*2,1);
            ylim = get(gca,'YLim')
            xlim('manual')
            xlimvec = [min(Xmatrix(:,1))-0.1,max(Xmatrix(:,1))+0.1];
            xlim(xlimvec);
            qscale = 1;
            hforce = quiver(Xmatrix(1,1),0,force(1),0,qscale);
            
            xlimvec = xlim;
            xvec = xlimvec(1):diff(xlimvec)/100:xlimvec(2);
            dBdx = -3 * xvec * coil.u0 * coil.I * coil.R^2./(2*(xvec.^2 + coil.R^2).^(5/2));
            Bscale = 10;
            plot(xvec,Bscale*dBdx,'-g')
            ylabel('dB/dx')
            set(gca,'YTickLabel',[])
            
            dt = diff(tvector);
            timetext = text(X(1),coil.R,'t = 0s');
            % loop through the simulation and update the visualization
            for tidx = 2:length(tvector)
                drawmagnet(Xmatrix(tidx,1),0,magnet.R*4,magnet.R*2,1,hmagnet);
                %disp(tvector(tidx))
                set(hforce,'XData',Xmatrix(tidx,1),'UData',force(tidx));
                timestr = strcat('t = ',num2str(tvector(tidx), '%3.1f'),' s');
                set(timetext,'String',timestr);
                
                % time delay for visualization
                pause(dt(tidx-1)/20)
            end
        end
        
    case 2
        disp('1 magnet 2 coils')
        coilarray(1) = coil;
        coilarray(2) = coil;
        coilarray(2).x = 1;
        coilarray(2).I = 50;
        
        x0_mag = .8;
        v0 = 0;
        X = [x0_mag v0];
        
        
        % Timing
        % ------
        tStart = 0;  % simulation start time (sec)
        tFinal = 100; % simulation start time (sec)
        
        
        graphicalsimualtion = true;
        %% Physics simulator
        % Set up the magnets starting config;
        Xinitial = X(end,:);
        
        tspan = [tStart tFinal];
        
        % solve for posn using ode45
        [tvector, Xmatrix] = ode45(@(t,state) compute_magnet_derivatives_2coils(t,state,magnet,coilarray), tspan, Xinitial);
        sol = ode45(@(t,state) compute_magnet_derivatives_2coils(t,state,magnet,coilarray), tspan, Xinitial);
        
        [xxdot,xdotxddot] = deval(sol,tvector);
        
        force = xdotxddot(2,:)/magnet.mass;
        
        figure
        subplot(3,1,1)
        plot(tvector,Xmatrix(:,1));
        title('Position')
        subplot(3,1,2)
        plot(tvector,Xmatrix(:,2));
        title('Velocity')
        
        subplot(3,1,3)
        plot(tvector,xdotxddot(2,:));
        title('Acceleration')
        xlabel('time (sec)')
        if graphicalsimualtion
            
            h1 = figure;
            hold on
            hcoil(1) = rectangle('Position',[coilarray(1).x-coilarray(1).R/8,-coilarray(1).R,coilarray(1).R/4,coilarray(1).R*2],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[coilarray(1).I/coilarray(1).I 0 0]);
            hcoil(2) = rectangle('Position',[coilarray(2).x-coilarray(2).R/8,-coilarray(2).R,coilarray(2).R/4,coilarray(2).R*2],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[coilarray(2).I/coilarray(2).I 0 0]);
            
            axis equal
            
            hmagnet = drawmagnet(Xmatrix(1,1),0,magnet.R*4,magnet.R*2,1);
            ylim = get(gca,'YLim')
            xlim('manual')
            
            %xlimvec = [min(Xmatrix(:,1))-0.1,max(Xmatrix(:,1))+0.1];
            %xlim(xlimvec);
            qscale = 1;
            hforce = quiver(Xmatrix(1,1),0,force(1),0,qscale);
            
            xlimvec = xlim;
            xvec = xlimvec(1):diff(xlimvec)/100:xlimvec(2);
            I = [coilarray(:).I];
            R = [coilarray(:).R];
            x = [coilarray(:).x];
            dBdx(1,:) = -3 * xvec * coilarray(1).u0 * I(1) .* R(1).^2./(2*((xvec-x(1)).^2 + R(1)^2).^(5/2));
            dBdx(2,:) = -3 * xvec * coilarray(2).u0 * I(2) .* R(2).^2./(2*((xvec-x(2)).^2 + R(2)^2).^(5/2));
            Bscale = 10;
            %plot(xvec,Bscale*sum(dBdx,1),'-g')
            ylabel('dB/dx')
            set(gca,'YTickLabel',[])
            
            dt = diff(tvector);
            timetext = text(X(1),coil.R,'t = 0s');
            % loop through the simulation and update the visualization
            for tidx = 2:length(tvector)
                drawmagnet(Xmatrix(tidx,1),0,magnet.R*4,magnet.R*2,1,hmagnet);
                %disp(tvector(tidx))
                set(hforce,'XData',Xmatrix(tidx,1),'UData',force(tidx));
                timestr = strcat('t = ',num2str(tvector(tidx), '%3.1f'),' s');
                set(timetext,'String',timestr);
                
                % time delay for visualization
                pause(dt(tidx-1)/20)
            end
        end
        
    case 3
        disp('1 magnet 2 coils, 1 variable current - trajectory following')
        
        % controller gains
        kd = 0.027;
        kp = 0.02; % 0.01 good if no kd
        
        coilarray(1) = coil;
        coilarray(2) = coil;
        coilarray(2).x = 1;
        coilarray(2).I = 50;
        
        % unpack variables for coil 2
        u0 = 4*pi*10^-7;
        R2 = coilarray(2).R;
        I2 = coilarray(2).I;
        L = coilarray(2).x;
        
        % set start and end positions
        xstart = 0.5; %start at equilibrium
        xgoal = 0.2;
        vzero = 0; 
        
        timetogoal = 2; %sec
        endtime = timetogoal*2; % how long the simulation will run for

        % feedback loop update frequency
        freq = 60;
        
        % generate trajectory
        [Q,QD,QDD] = jtraj(xstart,xgoal,freq*timetogoal,vzero,vzero);
        qtime = linspace(0,timetogoal,freq*timetogoal); % time vector for trajectory
        
        %% plot trajectory
%         hold all
%         plot(qtime,Q,'--r')
%         plot(qtime,QD,'--r')
%         plot(qtime,QDD,'--r')
%         legend({'Q';'QD';'QDD'})
%         title('Desired Trajectory')
        
        %% initialize state storage vectors
        tvector = 0;
        Xmatrix = [xstart,vzero];
        avector = 0;
        Ihistory = zeros(1,freq*endtime);
        Ihistory(1) = coilarray(1).I;
        
        for idx = 2:endtime*freq
            xnow = Xmatrix(end,1);
            vnow = Xmatrix(end,2);
            
            if idx > freq*timetogoal
                xdes = xgoal;
                vdes = 0;
            else
                xdes = Q(idx);
                vdes = QD(idx);
            end
            
            F1 = kd*(xdes-xnow) + kp*(vdes - vnow);
            I1 = ((F1 + u1*3*(xnow-L)*u0*I2*R2^2/(2*((xnow-L)^2 + R2^2)^(5/2)))*(2*(xnow^2+R1^2)))/(-3*xnow*u0*R1^2);
            Ihistory(idx) = I1;
            coilarray(1).I = I1;
            % run ode45 solve displacement for timestep
            tspan = [idx-1, idx]/freq;
            [tvectortemp, Xmatrixtemp] = ode45(@(t,state) compute_magnet_derivatives_2coils(t,state,magnet,coilarray), tspan, Xmatrix(end,:));
            sol = ode45(@(t,state) compute_magnet_derivatives_2coils(t,state,magnet,coilarray), tspan, Xmatrix(end,:));
        
            [xxdot,xdotxddot] = deval(sol,tvectortemp);
        
            forcetemp = xdotxddot(2,:)/magnet.mass;
        
            
            
            % store states
            tvector = [tvector;tvectortemp];
            Xmatrix = [Xmatrix;Xmatrixtemp];
            
            avector = [avector,xdotxddot(2,:)];
        end
        
        figure
        subplot(3,1,1)
        plot(tvector,Xmatrix(:,1));
        hold on
        plot(qtime,Q,'--r')
        title('Position')
        subplot(3,1,2)
        hold on
        plot(tvector,Xmatrix(:,2));
        plot(qtime,QD,'--r')
        title('Velocity')
        subplot(3,1,3)
        hold on
        plot(tvector,avector);
        plot(qtime,QDD,'--r')
        title('Acceleration')
        
        xlabel('time (sec)')
        
        %% plot current value
        figure
        itime = 0:1/freq:endtime;
        plot(itime(1:end-1),Ihistory)
        xlabel('time(sec)')
        ylabel('Current through coil 1 (A)')
        
        %% visualize
        if graphicalsimualtion
            
            h1 = figure;
            hold on
            hcoil(1) = rectangle('Position',[coilarray(1).x-coilarray(1).R/8,-coilarray(1).R,coilarray(1).R/4,coilarray(1).R*2],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[coilarray(1).I/coilarray(1).I 0 0]);
            hcoil(2) = rectangle('Position',[coilarray(2).x-coilarray(2).R/8,-coilarray(2).R,coilarray(2).R/4,coilarray(2).R*2],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[coilarray(2).I/coilarray(2).I 0 0]);
            
            axis equal
            
            hmagnet = drawmagnet(Xmatrix(1,1),0,magnet.R*4,magnet.R*2,1);
            ylim = get(gca,'YLim')
            xlim('manual')
            
            % plot start and goal
            %ylimvec = ylim;
            plot([xstart,xstart],ylim,'--g')
            plot([xgoal,xgoal],ylim,'--r')
            
            %xlimvec = [min(Xmatrix(:,1))-0.1,max(Xmatrix(:,1))+0.1];
            %xlim(xlimvec);
            qscale = 1;
            hforce = quiver(Xmatrix(1,1),0,avector(1),0,qscale);
            
            xlimvec = xlim;
            xvec = xlimvec(1):diff(xlimvec)/100:xlimvec(2);
            I = [coilarray(:).I];
            R = [coilarray(:).R];
            x = [coilarray(:).x];
            dBdx(1,:) = -3 * xvec * coilarray(1).u0 * I(1) .* R(1).^2./(2*((xvec-x(1)).^2 + R(1)^2).^(5/2));
            dBdx(2,:) = -3 * xvec * coilarray(2).u0 * I(2) .* R(2).^2./(2*((xvec-x(2)).^2 + R(2)^2).^(5/2));
            Bscale = 10;
            %plot(xvec,Bscale*sum(dBdx,1),'-g')
            ylabel('dB/dx')
            set(gca,'YTickLabel',[])
            
            dt = diff(tvector);
            timetext = text(Xmatrix(1,1),coil.R,'t = 0s');
            % loop through the simulation and update the visualization
            for tidx = 2:10:length(tvector)
                drawmagnet(Xmatrix(tidx,1),0,magnet.R*4,magnet.R*2,1,hmagnet);
                %disp(tvector(tidx))
                set(hforce,'XData',Xmatrix(tidx,1),'UData',avector(tidx));
                timestr = strcat('t = ',num2str(tvector(tidx), '%3.1f'),' s');
                set(timetext,'String',timestr);
                
                % time delay for visualization
                %pause(dt(tidx-1)/20)
            end
        end
        
        
    otherwise
        
        
        disp('Invalid simMode')
        
end