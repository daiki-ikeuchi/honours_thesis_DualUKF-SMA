## honours_thesis_DualUKF-SMA
# This repository contains the Matlab and Simulink files for Honours thesis completed at the University of Sydney.
# The thesis is available from https://www.researchgate.net/publication/336730812_Unscented_Kalman_Filter_based_Self-sensing_Control_of_Shape_Memory_Alloy_Actuator

Abstract:
A shape memory alloy actuator is widely used in various engineering fields due to its large force-to-weight ratio, 
large displacement, compact size and noiseless operation. However, the use of the actuator still remains uncommon 
compared with conventional actuators, especially as an active actuator. One of the reasons for such observation is attributed 
to the control difficulty due to the nonlinear and hysteresis nature of shape memory alloy and limited state measurement 
capabilities. Despite the great amount of research on the nonlinear and hysteresis sides, there has not been much research 
focus on the improvement of uncertainty associated with the lack of measured state variables. In this thesis, the author 
proposes to design a dual Unscented Kalman Filter model-based state and parameter estimator for a simple spring-biased SMA 
wire actuator and investigate the numerical feasibility of the estimation algorithms to overcome the aforementioned control 
issue. The findings of this study are expected to contribute to the shape memory alloy actuator control community to raise 
the awareness of significance of mitigating state and parameter uncertainties in terms of high-precision control through the 
use of estimation algorithms.

How to run:
1. Run Openloop_SMA_wire_param.slx (mimicking the actual experiments) or use the experimental data in your hand
2. Run main.m
3. Run Dual_UKF.slx
