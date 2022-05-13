function Vel = EncoderVelocity(position,time)
binWinTime = 2;
edgesPos = 0:1:abs(ceil(Velocity(end)));
VRposition =  discretize(position,edgesPos);
edgesTime = 0:binWinTime:time(end); %Bin width on the track
VRtimebin = discretize(time,edgesTime);
% Checks binning routine
checkBin = find(diff(VRtimebin)>1);
count = 0;
while checkBin>0
    for check = 1:size(checkBin,1)
        checkValue = checkBin(check)+1; %Adds 1 to the location of the skipped cell
        VRtimebin(checkValue,1) = VRtimebin(checkValue,1)-1;
    end
    checkBin = find(diff(VRtimebin)>1);
    count = count+1;
end
disp(['Binning routine adjusted: ' num2str(count)]);
clear count
% Checks for NaN values
nanValues = find(isnan(VRtimebin));
VRtimebin(isnan(VRtimebin)) = VRtimebin(nanValues(1)-1,1); % Assigns last value to NaN

dP = diff(position);
dP = vertcat(0,dP);
for i = 1:VRtimebin(end)
    findBinT = find(VRtimebin==i);
    sumPosition = sum(dP(findBinT,1));
    Vel(i,1) = i;
    Vel(i,2) = sumPosition;
end