# Free Energy Calculator
This repository includes the code to calculate the free energies and rate constants for an elementary surface reaction, adsorption desorption processes, or a simple vacuum phase reaction. 

Things to look out for
======================
* This code uses a wave number cutoff of 100/cm. Modify it if needed on FEGeneralProcedures.f90 file, line 90 and 91.
* Change the temperature range and interval you want to apply for in FEMain.f90 file, line 26 and 27.

Input File formats
==================
* For stationary point surface species, input files should be named as reac-1,2,3... or prod-1,2,3..etc. A sample format is shown below. Here the first line denotes the SCF energy in eV and the rest are associated wavenumbers (1/cm) for that species.

-421.24267709277507<br/>
168.823811<br/>
115.320187<br/>
19.695854<br/>
98.440556<br/>
125.522112<br/>
139.202853<br/>
197.017132<br/>
240.897136<br/>
343.524635<br/>
521.019461<br/>
684.776047<br/>
825.498616<br/>
880.318726<br/>
995.623088<br/>
1019.004489<br/>
1046.426720<br/>
1149.879318<br/>
1209.489890<br/>
1274.645908<br/>
1296.696516<br/>
1348.858457<br/>
1378.499895<br/>
1433.703007<br/>
1438.319513<br/>
2938.764063<br/>
2987.221790<br/>
3002.921773<br/>
3069.432205<br/>
3451.952704<br/>
3674.108309<br/>

* For a saddle point, the file name should be ts-1,2,3...etc. File format is same as before, do not include the imaginary wavenumber value in the input file.

* For a gas phase species, the file name should be Greac-1,2,3.. or Gprod-1,2,3... File format for a gas phase species is shown below. Here, the number on the first line represents the Zero Point Energy in eV, the sencod line represents the molecular weight of the species in kg/mole, and the following values in 2 columns represent temperature and total partition function (q) at the correspoinding temperature, respectively.  

-7.14516199804438<br/>
0.275693597480214<br/>
0.00202<br/>
273.00    1.3865793073029E+05<br/>
298.00    1.8842233198990E+05<br/>
323.00    2.4979530953306E+05<br/>
348.00    3.2426704082714E+05<br/>
373.00    4.1338658283972E+05<br/>
398.00    5.1875972043667E+05<br/>
423.00    6.4204704705429E+05<br/>
448.00    7.8496224976782E+05<br/>
473.00    9.4927058092918E+05<br/>
498.00    1.1367875061393E+06<br/>
523.00    1.3493775236479E+06<br/>
548.00    1.5889531535192E+06<br/>
573.00    1.8574740963278E+06<br/>
598.00    2.1569465611021E+06<br/>
623.00    2.4894227611095E+06<br/>
648.00    2.8570005742771E+06<br/>
673.00    3.2618233629226E+06<br/>
698.00    3.7060799453266E+06<br/>
723.00    4.1920047097311E+06<br/>
748.00    4.7218778597393E+06<br/>
773.00    5.2980257788972E+06<br/>
798.00    5.9228215014833E+06<br/>
823.00    6.5986852761999E+06<br/>
848.00    7.3280852095055E+06<br/>
873.00    8.1135379756987E+06<br/>
898.00    8.9576095814889E+06<br/>
923.00    9.8629161736075E+06<br/>
948.00    1.0832124878966E+07<br/>
973.00    1.1867954667889E+07<br/>
