function [] = Plot_Solution(solution, n)
    
    tt              = [];
	Altitude        = [];
	Longitude       = [];
	Latitude        = [];
	Speed           = [];
	Path_Angle      = [];
	Heading_Angle   = [];
	Mass            = [];

	Angle_Attack    = [];
	Bank_Angle      = [];
	Mass_Flow       = [];
    
    %Obtainting the interpolated solution from polynomials
	for i=1:n
        t               = linspace(solution(i).t0,solution(i).tf,solution(i).tf*100);
		tt              = [tt t];
		Altitude        = [Altitude; speval(solution(i),'X',1,t)];
		Longitude       = [Longitude; speval(solution(i),'X',2,t)];
		Latitude        = [Latitude; speval(solution(i),'X',3,t)];
		Speed           = [Speed; speval(solution(i),'X',4,t)];
		Path_Angle      = [Path_Angle; speval(solution(i),'X',5,t)];
		Heading_Angle   = [Heading_Angle; speval(solution(i),'X',6,t)];
		Mass            = [Mass; speval(solution(i),'X',7,t)];

		Angle_Attack    = [Angle_Attack; speval(solution(i),'U',1,t)];
		Bank_Angle      = [Bank_Angle; speval(solution(i),'U',2,t)];
		Mass_Flow       = [Mass_Flow; speval(solution(i),'U',3,t)];
    end

    figure1=figure;
    axes1 = axes('Parent',figure1,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    box(axes1,'on');
    grid(axes1,'on');
    hold(axes1,'all');
    % Create plot
    plot(tt,Altitude./1000,'LineWidth',3,'Color',[1 0 0],'Parent',axes1,...
        'DisplayName','Altitude');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes1,'show');

    figure
    geoplot(rad2deg(Latitude), rad2deg(Longitude), 'LineWidth',3,'Color',[1 0 0]);
    
    figure2=figure;
    axes2 = axes('Parent',figure2,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Geographic Coordinates (degrees)');
    box(axes2,'on');
    grid(axes2,'on');
    hold(axes2,'all');
    % Create plot
    plot(tt,rad2deg(Latitude),'LineWidth',3,'Color',[1 0 0],'Parent',axes2,...
        'DisplayName','Latitude');
    plot(tt,rad2deg(Longitude),'LineWidth',3,'Color',[0 1 0],'Parent',axes2,...
        'DisplayName','Longitude');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes2,'show');


    figure3=figure;
    axes3 = axes('Parent',figure3,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Speed (km/s)');
    box(axes3,'on');
    grid(axes3,'on');
    hold(axes3,'all');
    % Create plot
    plot(tt,Speed,'LineWidth',3,'Color',[1 0 0],'Parent',axes3,...
        'DisplayName','Speed');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes3,'show');

    figure4=figure;
    axes4 = axes('Parent',figure4,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Angles (degrees)');
    box(axes4,'on');
    grid(axes4,'on');
    hold(axes4,'all');
    % Create plot
    plot(tt,rad2deg(Path_Angle),'LineWidth',3,'Color',[1 0 0],'Parent',axes4,...
        'DisplayName','Path Angle');
    plot(tt,rad2deg(Heading_Angle),'LineWidth',3,'Color',[0 1 0],'Parent',axes4,...
        'DisplayName','Heading Angle');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes4,'show');

    figure5=figure;
    axes5 = axes('Parent',figure5,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Angles (degrees)');
    box(axes5,'on');
    grid(axes5,'on');
    hold(axes5,'all');
    % Create plot
    plot(tt,rad2deg(Angle_Attack),'LineWidth',3,'Color',[1 0 0],'Parent',axes5,...
        'DisplayName','Angle of Attack');
    plot(tt,rad2deg(Bank_Angle),'LineWidth',3,'Color',[0 1 0],'Parent',axes5,...
        'DisplayName','Bank Angle');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes5,'show');


    figure6=figure;
    axes6 = axes('Parent',figure6,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Mass (kg)');
    box(axes6,'on');
    grid(axes6,'on');
    hold(axes6,'all');
    % Create plot
    plot(tt,Mass,'LineWidth',3,'Color',[1 0 0],'Parent',axes6,...
        'DisplayName','Mass');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes6,'show');

    figure7=figure;
    axes7 = axes('Parent',figure7,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Mass Flow Rate (kg/s)');
    box(axes7,'on');
    grid(axes7,'on');
    hold(axes7,'all');
    % Create plot
    plot(tt,Mass_Flow,'LineWidth',3,'Color',[1 0 0],'Parent',axes7,...
        'DisplayName','Fuel Injection Rate');
    % Create title
    title('Minimum Time to Climb (h-Method)','FontSize',20);
    % Create legend
    legend(axes7,'show');

end

