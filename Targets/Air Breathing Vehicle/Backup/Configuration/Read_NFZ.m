function [Alt2,Lon2,Lat2,Radius,n_nfz,error] = Read_NFZ(config)
    xDoc = xmlread(config);                       
    allRadius = getElementsByTagName(xDoc,'Radius');
    allOrigin = getElementsByTagName(xDoc,'Origin');
	error = 0;
	Alt2 = [];
	Lon2 = [];
	Lat2 = [];
	Radius = [];
	
    %Total number of No Fly Zones
	n_nfz = allRadius.getLength;
    total_nfz = n_nfz-1;

    %Extract and store the details of the No Fly Zones
    for i=0:total_nfz
        Radius_Node = item(allRadius,i);
        Radius(i+1) = str2double(getFirstChild(Radius_Node).getData);
        error = error + Check_Limit(Radius(i+1), 200, 0);

        Origin_Node = item(allOrigin, i);
        Child_Origin = getNextSibling(getFirstChild(Origin_Node));
        Lat2(i+1) = str2double(Child_Origin.getFirstChild.getData);
        Child_Origin = getNextSibling(getNextSibling(Child_Origin));
        Lon2(i+1) = str2double(Child_Origin.getFirstChild.getData);
        Child_Origin = getNextSibling(getNextSibling(Child_Origin));
        Alt2(i+1) = str2double(Child_Origin.getFirstChild.getData);
        error = error + Check_Limit(Alt2(i+1), 80e3, 200);
    end
end

