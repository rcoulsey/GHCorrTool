**This code is incomplete. Spatial correction should be added. 
For now, it is an unconstrained correlation of the windowed time series

Class Sensor:
variable names:
Z_30 Vertical component time series up to 30Hz
Power_30 Power spectra up to 30Hz for sensor
KR Spacial correction coefficients for Coherance Calculations

The code has sensor object constructors for 9 stations written in, change this for your array.

For each ring, coherance is calculated by the sensor objects you pass to it. 
Currently set up for 3 rings with 3, 6, and 6 sensors as the rings move out:
change for use case. 
