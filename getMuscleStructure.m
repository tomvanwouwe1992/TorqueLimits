function [MusclesOI,CoordinatesOI] = getMuscleStructure( Model_OS,CoordinatesOI_Names )
% Print out the joints of the model
Nr_Coordinates = Model_OS.getCoordinateSet.getSize;
for i = 1:Nr_Coordinates
    CoordinateNames(i).name = Model_OS.getCoordinateSet.get(i-1).getName.toCharArray';
end

% We are interested in muscles that actuate HIP(3) - KNEE - ANKLE of the
% right leg
for i = 1:length(CoordinatesOI_Names)
    for j = 1:length(CoordinateNames)
        if CoordinateNames(j).name == string(CoordinatesOI_Names(i));
            CoordinatesOI(i).index = j;
            
            break
        end
    end
end

% % Print out the coördinates in the Joints Of Interest
% for i = 1:length(CoordinatesOI)
%     index = CoordinatesOI(i).index;
%     sizeCoordinateSet = Model_OS.getJointSet.get(index-1).getCoordinateSet.getSize;
%     for j = 1:sizeCoordinateSet
%         JointsOI(i).CoordinateSet(j).index = j;
%         JointsOI(i).CoordinateSet(j).name = Model_OS.getJointSet.get(index-1).getCoordinateSet.get(j-1).getName;
%     end
% end

state = Model_OS.initSystem;
Model_OS.equilibrateMuscles(state);
Nr_Muscles = Model_OS.getMuscles.getSize;
incM = 0;
MuscleNr = 1;
for m = 1:Nr_Muscles
    muscle = Model_OS.getMuscles.get(m-1);
    Coordinate_Nr = 1;
    for i = 1:length(CoordinatesOI)
        coordinateIndex = CoordinatesOI(i).index;
        coordinate = Model_OS.getCoordinateSet.get(coordinateIndex-1);
        if abs(muscle.computeMomentArm(state, coordinate)) > 0.0001
            MusclesOI(MuscleNr).name = muscle.getName;
            MusclesOI(MuscleNr).index = m;
            MusclesOI(MuscleNr).coordinate(Coordinate_Nr).coordinateName = coordinate.getName.toCharArray';
            MusclesOI(MuscleNr).coordinate(Coordinate_Nr).coordinateIndex = coordinateIndex;
            incM = 1;
            Coordinate_Nr = Coordinate_Nr + 1;
        end
    end
    if incM == 1
        MuscleNr = MuscleNr + 1;
    end
    inc = 0;
end



% Put the Muscle Properties into structure
MusclesOI = getMuscleProperties(Model_OS,MusclesOI);



end

