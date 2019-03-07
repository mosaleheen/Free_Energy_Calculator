# Free_Energy_Calculator
This repository includes the code to calculate the free energies and rate constants for an elementary surface reaction, adsorption desorption processes, or a simple vacuum phase reaction. 

Things to look out for
======================
1. Author does not guarantee that this package is free from error. Neither does he accept responsibility
   for any loss or damage that results from its use.
1. This code uses a wave number cutoff of 100/cm. Modify it if needed on FEGeneralProcedures.f90 file, line 90 and 91.
2. Chage the temperature and intervals you want to apply for in FEMain.f90 file, line 26 and 27.

Input File formats
==================
1. For stationary point surface species, input files can be named as reac-1,2,3... or prod-1,2,3...depending on what they are. Sample format is shown below. Here the first line denotes the SCF energy in eV and the rest are associated wavenumbers (1/cm) for that species.

-421.24267709277507
168.823811
115.320187
19.695854
98.440556
125.522112
139.202853
197.017132
240.897136
343.524635
521.019461
684.776047
825.498616
880.318726
995.623088
1019.004489
1046.426720
1149.879318
1209.489890
1274.645908
1296.696516
1348.858457
1378.499895
1433.703007
1438.319513
2938.764063
2987.221790
3002.921773
3069.432205
3451.952704
3674.108309

2. For a saddle point, the file name should be ts-1,2,3...etc. File format is same as before, except do not include the imaginary wavenumber value in the input file.
3. For a gas phase species, the file name should be Greac-1,2,3.. or Gprod-1,2,3... File format for a gas phase species is shown below. Here, the number on the first line represents the Zero Point Energy in eV, the sencod line represents the molecular weight of the species in kg/mole, and the following values in 2 columns represents temperature and partition function (q) at the that temperature.  

-7.14516199804438
0.275693597480214
0.00202
273.00    1.3865793073029E+05
298.00    1.8842233198990E+05
323.00    2.4979530953306E+05
348.00    3.2426704082714E+05
373.00    4.1338658283972E+05
398.00    5.1875972043667E+05
423.00    6.4204704705429E+05
448.00    7.8496224976782E+05
473.00    9.4927058092918E+05
498.00    1.1367875061393E+06
523.00    1.3493775236479E+06
548.00    1.5889531535192E+06
573.00    1.8574740963278E+06
598.00    2.1569465611021E+06
623.00    2.4894227611095E+06
648.00    2.8570005742771E+06
673.00    3.2618233629226E+06
698.00    3.7060799453266E+06
723.00    4.1920047097311E+06
748.00    4.7218778597393E+06
773.00    5.2980257788972E+06
798.00    5.9228215014833E+06
823.00    6.5986852761999E+06
848.00    7.3280852095055E+06
873.00    8.1135379756987E+06
898.00    8.9576095814889E+06
923.00    9.8629161736075E+06
948.00    1.0832124878966E+07
973.00    1.1867954667889E+07
