0,0.1,0.0,0.2         !# computational domain [XA,XB][YA,YB]
200,411        !# number of nodes !even number1
2,30					                  !# terminal time,  number_write
3,1,1,2                                  !# boundry condition 1=noslip 2=slip,3=axisymmetric,4=zero gradient,LEF,RIG,BOT,TOP
996,1.0e-3                                !# water:density  viscosity 1.0087e-3  
1.225,1.77E-5                             !# air  :density  viscosity 1.81e-5 
1e-3                                      !# mass residual
0.073                                     !# epiron
90,90                            	      !# thetaA,thetaR
0.1,1.35e-3,7                             !# h0,r0,deg
4.e-4	                                    !#mindt
1.8,15                                   ！amplitude,frequency
3                                          !NS_correction iteration max cycle number
2                                          !set_compute

  




