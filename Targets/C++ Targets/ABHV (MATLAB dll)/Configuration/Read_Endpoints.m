function [Alt0,Lon0,Lat0,Altf,Lonf,Latf,err] = Read_Endpoints(config)
    try
       xDoc = xmlread(config); 
    catch
       err = -1;
       error('Failed to read XML file %s.', config);
    end
                          
    Takeoff_Position = getElementsByTagName(xDoc,'Takeoff_Position');
    Landing_Position = getElementsByTagName(xDoc,'Landing_Position');
    err = 0;
    
    %Extract Initial (Takeoff) Position of the Vehicle
    Coordinates_Takeoff = item(Takeoff_Position, 0);
    Child_Takeoff = getNextSibling(getFirstChild(Coordinates_Takeoff));
    Lat0 = str2double(Child_Takeoff.getFirstChild.getData);
    Child_Takeoff = getNextSibling(getNextSibling(Child_Takeoff));
    Lon0 = str2double(Child_Takeoff.getFirstChild.getData);
    Child_Takeoff = getNextSibling(getNextSibling(Child_Takeoff));
    Alt0 = str2double(Child_Takeoff.getFirstChild.getData);
    err = err + Check_Limit(Alt0, 80e3, 80);

    %Extract Final (Landing) Position of the Vehicle
    Coordinates_Landing = item(Landing_Position, 0);
    Child_Landing = getNextSibling(getFirstChild(Coordinates_Landing));
    Latf = str2double(Child_Landing.getFirstChild.getData);
    Child_Landing = getNextSibling(getNextSibling(Child_Landing));
    Lonf = str2double(Child_Landing.getFirstChild.getData);
    Child_Landing = getNextSibling(getNextSibling(Child_Landing));
    Altf = str2double(Child_Landing.getFirstChild.getData);
    err = err + Check_Limit(Altf, 80e3, 80);
end

