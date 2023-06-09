function [Alt1,Lon1,Lat1,Speed,n_wp,error] = Read_Waypoints(config)
    xDoc = xmlread(config);                       
    allSpeed = getElementsByTagName(xDoc,'Speed');
    allCoordinates = getElementsByTagName(xDoc,'Coordinates');
    error = 0;
	Alt1 = [];
	Lon1 = [];
	Lat1 = [];
	Speed = [];
	
    %Total number of Waypoints
	n_wp = allSpeed.getLength;
    total_waypoints = n_wp-1;

    %Extract and store the details of the Waypoints
    for i=0:total_waypoints
        Speed_Node = item(allSpeed,i);
        Speed(i+1) = str2double(getFirstChild(Speed_Node).getData);
        error = error + Check_Limit(Speed(i+1), 2800, 120);

        Coordinates_Node = item(allCoordinates, i);
        Child_Coordinates = getNextSibling(getFirstChild(Coordinates_Node));
        Lat1(i+1) = str2double(Child_Coordinates.getFirstChild.getData);
        Child_Coordinates = getNextSibling(getNextSibling(Child_Coordinates));
        Lon1(i+1) = str2double(Child_Coordinates.getFirstChild.getData);
        Child_Coordinates = getNextSibling(getNextSibling(Child_Coordinates));
        Alt1(i+1) = str2double(Child_Coordinates.getFirstChild.getData);
        error = error + Check_Limit(Alt1(i+1), 80e3, 80);
    end
end

