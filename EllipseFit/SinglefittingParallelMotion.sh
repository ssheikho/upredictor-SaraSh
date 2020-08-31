#!/bin/bash  
# int-or-string.sh

make

declare -i reachCount
# The script will treat subsequent occurrences of "number" as an integer.	
	

echo "ID, ReachSample #, Phase, n, Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, inPtSx, inPtSy, inPtSz, inPtTx, inPtTy, inPtTz, inEF start Theta, inEF target Theta, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Ellipse2D CX, Ellipse2D CY, Ellipse2D A, Ellipse2D B, Ellipse2D B/A, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" > ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv

#echo "ID, Reach Count, Phase, size (n), snippet %, inPts(1)/predPts(2), 3dPtsX, 3dPtsY, 3dPtsZ, fitErrResidualsN train(1)/test(2) " >  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/FittedEllipse_inVsOut_PathPts-P2.txt	


reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_151_344_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_151_392_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_515_710_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_515_717_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1042_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1080_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1211_1410_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1211_1434_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_2164_2342_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_2463_2654_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4598_4797_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4900_5085_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_5200_5392_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_5810_5999_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_6104_6298_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_6396_6559_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_8235_8437_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_8235_8464_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_8588_8787_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_8588_8807_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_8924_9122_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_9247_9434_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_9848_10035_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_10156_10350_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_10443_10631_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_12183_12366_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_12183_12420_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_12520_12725_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_12819_13017_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_13162_13335_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_13716_13894_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_14004_14169_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_15853_16175_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_16265_16486_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_16584_16783_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_16892_17059_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_17168_17337_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_17168_17374_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_17458_17650_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_17768_17942_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_18033_18208_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_19751_20003_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_20127_20324_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_20127_20334_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_20472_20671_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_20472_20710_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_20812_21005_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_21104_21304_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_21104_21307_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_21384_21567_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_21674_21839_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_23596_23787_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24237_24427_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_25122_25303_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_25420_25604_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_27509_27818_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_27927_28109_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_28531_28722_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_28822_29288_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_29121_29288_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_29412_29590_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_31685_31903_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_32013_32208_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_32326_32505_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_32615_32779_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_32881_33056_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_33180_33332_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_33451_33635_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_35222_35486_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_35932_36123_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_36235_36414_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_36843_37031_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_37190_37376_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_39428_39627_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_39428_39631_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_39815_40015_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_39815_40025_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_40174_40395_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_40570_40760_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_40570_40768_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_41122_41305_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_41431_41627_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_43594_43808_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_44261_44457_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_44568_44743_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_44839_45008_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_45122_45310_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_45417_45581_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_45687_45849_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))

#echo "Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" >  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/JointFitting_Plane_PCAPts-P3.txt

reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_479_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1352_1552_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1352_1603_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1741_1988_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_4756_5063_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_5875_6107_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_6233_6471_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_6999_7258_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_10423_10681_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_10814_11082_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_11563_11786_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_13990_14265_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_15188_15468_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_15594_15831_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_16399_16678_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_20281_20549_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_20664_20968_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_21478_21760_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_24055_24387_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_24486_24824_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_26674_26977_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_30877_31174_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_31308_31602_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_32188_32456_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_36431_36712_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_37772_38061_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_40490_40838_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_41806_42122_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_42264_42584_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_281_594_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_1541_1841_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_1985_2308_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_2931_3195_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_7701_7953_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_8547_8810_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_13972_14284_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_18054_18334_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_18459_18732_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_19312_19609_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_24016_24303_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_24863_25148_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_30312_30617_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_33217_33636_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_35124_35434_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_36069_36382_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_2/2016_11_24_study_4_Part3_2_41682_42001_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))

#echo "Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" >  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/JointFitting_Plane_PCAPts-P10.txt

reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_834_1101_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_1250_1517_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_1705_1990_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2604_2859_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_3054_3318_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_3449_3719_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_3873_4136_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6474_6725_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_7251_7472_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_7607_7808_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_7971_8188_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_8373_8617_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_9051_9246_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_10489_10871_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_11240_11511_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_11645_11914_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_12073_12321_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_12431_12693_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_12797_12987_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_13160_13398_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_13837_14053_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_15956_16175_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_16296_16546_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_16675_16901_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_17029_17199_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_17366_17606_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_17736_17988_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_18089_18333_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_18439_18669_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_21445_21708_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_21854_22127_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_22236_22531_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_22611_22850_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_23003_23252_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_23360_23595_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_23720_23955_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_25659_25906_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_26026_26268_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_26392_26608_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_26820_27098_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_27196_27461_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_27572_27786_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_27981_28228_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_28329_28608_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_28708_28986_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_29381_29625_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_31252_31514_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_31632_31942_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_32064_32326_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_32446_32662_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_32827_33040_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_33222_33489_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_33595_33837_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_33963_34265_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_36294_36593_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_36723_37019_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_37156_37451_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_37564_37912_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_37998_38291_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_38444_38732_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_38836_39094_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_39234_39522_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_41619_41921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_42000_42283_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_42415_42693_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_43282_43562_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_43693_44002_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_44085_44484_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_44621_44910_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_46720_47009_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_47110_47364_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_47475_47656_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_47858_48125_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_48217_48457_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_48616_48824_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_49004_49278_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_49381_49668_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_49773_50057_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_293_590_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_728_1042_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_1131_1360_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_1578_1879_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_2025_2315_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_2443_2709_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_2844_3136_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_3247_3495_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_3646_3942_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_5337_5936_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_6074_6391_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_6512_6788_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_6935_7218_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_7377_7666_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_7802_8079_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_8701_8957_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_9081_9351_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_11745_12089_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_12178_12410_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_12626_12932_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_13048_13415_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_13529_13825_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_13964_14280_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_14381_14664_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_14802_15093_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_17019_17292_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_17434_17746_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_17879_18159_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_18324_18644_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_18769_19049_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_19200_19482_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_19644_19953_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_20047_20314_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_20444_20732_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_22565_22883_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_23006_23347_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_23431_23695_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_23884_24172_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_24281_24518_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_24700_24976_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_25147_25457_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_25573_25865_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_25981_26252_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_28347_28638_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_28752_29119_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_29209_29470_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_29671_29993_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_30101_30428_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_30537_30817_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_30970_31266_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_31394_31693_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_31805_32119_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_34077_34413_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_34503_34867_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_34942_35262_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_35395_35678_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_35791_36034_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_36204_36495_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_36655_36931_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_37054_37351_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_37467_37751_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_39643_39921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_40034_40376_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_40479_40688_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_40872_41128_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_41244_41552_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_41629_41931_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_42032_42289_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_42806_43102_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_44956_45232_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_45350_45670_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_45772_46001_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_46177_46461_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_46554_46868_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_46954_47215_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_47379_47640_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_47760_48022_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_2/2016_12_14_study_4_Part10_2_48136_48410_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))

#echo "Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv

reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_495_761_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_895_1077_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_1241_1433_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_1638_1878_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_2001_2213_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_2369_2696_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_2831_3063_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_3193_3408_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_3550_3760_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_5747_5998_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_6111_6286_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_6482_6670_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_6850_7091_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_7194_7460_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_7543_7773_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_7912_8147_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_8274_8522_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_8640_8846_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_10470_10680_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_10791_11000_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_11142_11314_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_11488_11711_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_11816_12005_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_12143_12362_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_12473_12680_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_12797_12978_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_13114_13318_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_15464_15645_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_15805_15972_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_16149_16376_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_16493_16682_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_16826_17067_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_17220_17470_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_17597_17807_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_17921_18148_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_20274_20486_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_20629_20804_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_20985_21200_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_21307_21513_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_21635_21830_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_21999_22197_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_22309_22527_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_22632_22882_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_24794_24947_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_25123_25283_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_25466_25660_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_25765_25942_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_26084_26310_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_26404_26582_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_26680_26882_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_26995_27200_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_29231_29460_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_29569_29721_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_29908_30119_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_30243_30480_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_30914_31115_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_31261_31466_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_31598_31809_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_33610_33825_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_34252_34451_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_34559_34796_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_35223_35421_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_35522_35743_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_35861_36086_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_37555_37921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_38049_38331_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_38404_38616_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_38868_39089_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_39210_39395_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_39570_39816_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_39942_40169_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_40302_40485_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_40672_40872_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_43108_43361_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_43794_44011_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_44116_44316_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_44752_44956_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_45372_45585_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_47305_47621_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_47687_47890_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_48021_48226_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_48343_48612_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_48784_48970_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_49078_49344_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_49488_49696_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_49878_50099_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_50244_50439_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_50562_50889_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_1091_1289_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_1447_1616_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_1817_2051_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_2142_2322_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_2489_2744_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_2845_3041_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_3191_3401_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_3528_3744_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_5189_5521_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_5650_5811_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_5955_6119_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_6301_6535_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_6654_6910_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_7317_7519_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_7655_7862_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_7974_8193_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_10026_10190_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_10339_10601_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_10683_10900_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_11005_11245_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_11330_11579_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_11698_11921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_12055_12297_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_12418_12641_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_14164_14454_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_14577_14746_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_14912_15193_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_15283_15487_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_15623_15880_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_16310_16542_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_16708_16920_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_17055_17278_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_18571_18926_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_19050_19319_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_19427_19734_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_20145_20402_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_20838_21042_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_21188_21381_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_21509_21717_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_23578_23744_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_23906_24050_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_24218_24417_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_24504_24711_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_25114_25300_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_25401_25583_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_25698_25886_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_27116_27441_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_27556_27838_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_27926_28187_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_28592_28818_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_29239_29449_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_29568_29771_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_29914_30112_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_32147_32368_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_32492_32666_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_32849_33080_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_33195_33453_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_33549_33804_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_33908_34154_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_34292_34517_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_34635_34831_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_36215_36548_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_36663_36921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_37019_37193_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_37392_37583_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_37700_37885_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_38027_38255_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_38365_38574_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_38689_38904_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_40411_40724_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_40830_41091_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_41483_41681_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_41801_42026_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_42437_42629_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_42734_42936_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_44672_45072_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_45310_45482_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_45652_45885_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_46037_46263_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_46367_46514_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_46652_46913_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_47022_47231_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_47343_47557_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_2/2016_12_14_study_4_Part12_2_47691_47980_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))

#echo "Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" >  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/JointFitting_Plane_PCAPts-P13.txt

reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_301_550_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_673_899_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_1028_1254_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_1394_1648_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_1772_2015_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_2131_2339_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_2487_2733_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_2872_3104_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_3229_3543_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_5074_5296_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_5401_5610_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_5746_5962_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_6091_6306_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_6419_6637_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_7137_7334_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_7451_7657_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_7771_7977_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_9531_9773_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_9872_10099_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_10223_10398_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_10589_10804_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_10928_11162_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_11595_11805_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_11923_12128_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_12237_12430_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_13992_14281_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_14368_14557_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_14751_14934_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_15172_15389_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_15548_15786_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_16248_16491_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_16621_16836_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_16958_17194_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_18470_18808_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_19089_19248_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_19399_19604_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_19751_19991_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_20369_20572_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_20713_20906_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_21034_21242_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_21356_21579_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_23600_23875_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_23934_24100_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_24328_24522_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_24712_24985_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_25093_25339_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_25467_25691_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_25937_26207_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_26337_26567_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_26690_26914_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_27991_28393_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_28569_28716_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_28869_29023_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_29303_29523_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_29939_30105_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_30502_30679_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_30782_30972_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_32448_32656_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_32754_33209_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_33366_33556_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_33976_34167_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_34543_34726_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_34844_35025_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_37054_37249_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_37357_37562_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_37664_37842_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_38043_38237_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_38327_38879_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_39025_39237_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_39345_39542_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_39652_39849_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_42015_42241_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_42356_42550_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_42719_42880_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_43037_43295_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_43630_43827_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_44228_44405_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_1/2016_12_14_study_4_Part13_1_44501_44685_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_229_499_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_863_1019_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_1175_1407_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_1760_1937_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_2350_2525_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_2637_2799_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_4863_5099_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_5483_5641_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_5794_6000_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_6395_6569_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_6667_6878_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_7291_7486_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_9517_9757_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_10176_10329_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_10744_10905_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_11315_11488_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_11567_11753_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_11863_12063_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_13888_14146_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_14487_14628_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_15045_15242_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_15634_15817_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_15894_16091_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_16159_16349_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_18055_18290_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_18345_18487_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_18915_19125_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_19462_19633_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_19724_19912_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_20256_20423_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_22114_22289_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_22394_22536_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_22976_23150_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_23275_23674_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_23817_23991_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_24084_24244_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_24350_24516_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_25546_25787_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_25883_26127_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_26482_26637_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_26813_26982_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_27152_27296_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_27434_27629_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_27846_28033_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_28167_28380_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_28509_28704_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_30490_30662_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_30807_30951_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_31107_31283_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_31455_31677_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_31988_32156_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_32269_32668_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_33932_34301_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_34497_34636_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_35070_35274_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_35345_35740_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_36145_36544_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_37586_37970_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_38351_38490_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_38871_39028_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_39426_39826_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_39922_40090_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_41214_41613_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_41783_41970_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_42116_42248_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_42395_42794_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_43179_43357_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_43447_43625_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_45001_45385_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_45764_45917_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_46360_46522_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_46634_46798_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_47143_47319_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part13_2/2016_12_14_study_4_Part13_2_47418_47590_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))

#echo "Eigenvalue_1, Eigenvalue_2, Eigenvalue3, v3_x, v3_y, v3_z, planeNToZ, planeYawZ, planeDipX, planeTiltY, Ellipse2D Alpha, Max AngularErrS, Min AngularErrs, HypErrShellAxesX, HypErrShellAxesY, HypErrShellAxesZ, fixedLEcenterOffset, fixedLEllipsoidAxesX, fixedLEllipsoidAxesY, fixedLEllipsoidAxesZ" >  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/JointFitting_Plane_PCAPts-P9.txt

reachCount=1

./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_362_566_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_705_899_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_1021_1194_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_1411_1621_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_2026_2202_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_2364_2549_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_2667_2822_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_3028_3212_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_5417_5621_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_5743_5880_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_6036_6195_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_6371_6579_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_7017_7184_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_7609_7790_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_7880_8057_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_9785_9995_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_10101_10269_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_10754_10932_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_11371_11595_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_12381_12564_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_14260_14453_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_14884_15026_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_15189_15369_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_15774_15960_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_18573_18682_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_18896_19041_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_19183_19338_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_19499_19678_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_19783_19988_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_20389_20548_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_20927_21083_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_22831_23046_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_23163_23309_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_23453_23599_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_24065_24262_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_24395_24576_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_24699_25098_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27332_27477_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27608_27749_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_28150_28321_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_28456_28648_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_28752_28921_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_29009_29408_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_31013_31258_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_32198_32331_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_32752_32915_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_33012_33411_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_34724_34925_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_35314_35469_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_35892_36039_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_36688_36855_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_36954_37119_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_39104_39247_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_39383_39549_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_39701_39874_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_39972_40115_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_41094_41258_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_43541_43691_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  
./EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_45160_45334_clean_phase1.csv $reachCount >> ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/Nominal\ Plane\ Fit\ in\ PCA/NominalPlaneFitPars.csv
 ((reachCount = reachCount +1))
  >>  ~/Dropbox/Reaching\ Study/Phase1/Experimental\ Data/EllipseFit/NewData/JointFitting_Plane_PCAPts-P9.txt	

           
