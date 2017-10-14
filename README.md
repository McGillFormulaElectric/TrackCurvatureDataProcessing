# TrackCurvatureDataProcessing
Track Curvature Data Processing Tool


1) create trackProcessingTool object
t = trackProcessingTool ;

2) add coordinates for track to be created to the object
t.X = xPoints ;
t.Y = yPoints ;

3) add event for which the track is used

4) optionally add units in which coordinates are given and if number of laps for endurance track

5) optionally specify step size used to generate radii and step size for output track

6) generate track
track = t.createTrack ;

7) a figure will be generated with calculated instantaneous radius and smoothed radius. If results 
are satisfactory save the track, otherwise try using a different filter parameter
