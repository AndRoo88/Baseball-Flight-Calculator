# Baseball-Flight-Calculator
Open Source Baseball Flight Simulator using RK4 numerical integration method. This project is designed to be modifiable and improvable see: https://baseballaero.com/UMBA

Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

Prerequisites
This is a python code. you will need python installed as well as the numpy and matplotlib.pyplot libraries installed I have found https://www.codecademy.com/articles/install-python as a useful link. Additionally, you can run the code from an ide like spyder. Spyder (via anaconda) has the libraries pre-installed.https://www.anaconda.com/distribution/

Once python and the requisite libraries are installed you can run the code.

insert the following ball variables in order with spaces between them #( for now)
x(in ft) y(in ft) z(in ft) vTot(mph) Theta(deg) Psi(deg) SpinRate(rpm) tilt(deg) gyro(deg)

x is the location left of the pitcher
y is the position from the rubber
z is the height of the ball
vTot is the ball speed
Theta is the ball's launch angle above horizontal
Psi is the heading angle in degrees right of pitcher
SpinRate is spin rate of ball
Tilt is the degrees of tilt on the clock face (i.e. 0 deg = 12:00 Tilt, 90 deg = 3:00, -90 deg = 9:00 Tilt, etc)
Gyro is the degrees away from the y direction

Two additional angles are asked for,. They are the angles from the spin axis to the ball logo.(not yet used)

This is built to enable pitchers, pitching coaches, or analysts gain a better understanding of ball flight behaviour and is built to be ever more accurate as better models and results become available

Versioning
We use SemVer for versioning. For the versions available, see the tags on this repository.

Authors
Andrew Smith - Seam Effect modeling and Initial work - USU EFDL
Michael Ressler - GitHub expertise -Diamond Kinetics(?)

Acknowledgments
Thanks to Barton Smith for enabling this work to get started, Alan Nathan for access to the excel file and answers to question and everyone else
etc
