# dynaPorousFoam

- c++ test folder is the place that I used to save the mid-term codes
- dynaPorousFoam
    This solver inherit the porous media method, where the input variables include D and F which are used in Dacy-Forschheimer equation, porosity which is used to indicate the position of net in the fluid domain. 
- dynaPorousFoamv2 
    This solver is developed for coupling with external FE solver. Code_Aster is my first option, of course, other FE solvers also have the potential to coupling. 
    The input variables include: (1) positions of each node on the netting, (2) surface element of net panel and (3) hydrodyanmic forces on each net panel.
    
