#!/bin/bash  
# int-or-string.sh
#make clean 
make

declare -i reachCount
			
# (A).Maya 
#		(A.1).inPts vicon

#P3_1	

#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_30877_31174_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_30877_31174_clean_phase1_ViconMrks.txt

#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_30877_31174_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_30877_31174_clean_phase2_ViconMrks.txt

#        ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_20664_20968_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_20664_20968_clean_phase1_ViconMrks.txt  

 #       ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_20664_20968_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_20664_20968_clean_phase2_ViconMrks.txt  


#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase1_ViconMrks.txt  

#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase2_ViconMrks.txt  


#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1352_1552_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase1_ViconMrks.txt 

#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_ViconMrks.txt 
        
        
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_37772_38061_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_37772_38061_clean_phase1_ViconMrks.txt

#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_37772_38061_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_37772_38061_clean_phase2_ViconMrks.txt


#		(A.2).inPts Ellipse-Fit

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_30877_31174_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_30877_31174_clean_phase1_Ef.txt 

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_30877_31174_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_30877_31174_clean_phase2_Ef.txt 

 #./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_20664_20968_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_20664_20968_clean_phase1_Ef.txt   

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_20664_20968_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_20664_20968_clean_phase2_Ef.txt   


#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt   

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt   


#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_1352_1552_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase1_Ef.txt  

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_Ef.txt  
        
        
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_37772_38061_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_37772_38061_clean_phase1_Ef.txt 

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_37772_38061_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_37772_38061_clean_phase2_Ef.txt 







# (A).Maya 
#		(A.1).inPts vicon
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_ViconMrks.txt

		
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_Ef_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_1352_1552_clean_phase2_Ef.txt  









 #       ./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef.txt  
  #      ./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_stusdy_4_Part3_1_2510_2749_clean_phase2_Ef.txt  





# (A). Fit elliptical model to inPts
make ../EllipseFit
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase1.txt

#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase2.txt





# (B).Maya - OutPTs
#./PathPlanning  outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase1.txt
#./PathPlanning  outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt

# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_253_453_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_453_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_253_453_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_253_453_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_253_453_clean_phase2_Ef.txt










#P3_1		example 2: 2510_2749	
# (A). Fit elliptical model to inPts
make ../EllipseFit
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.txt

#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2.txt





# (B).Maya - OutPTs
#./PathPlanning  outEfPts/2016_11_24_study_4_Part3_1_ 2510_2749_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.txt
#./PathPlanning  outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_ViconMrks.txt



#		(A.2).inPts Ellipse-Fit

#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef.txt

# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_2510_2749_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part3_1_2510_2749_clean_phase2_Ef.txt







#  Part9_- Example 1. 38812_39041

# (0). Fit elliptical model to inPts
make ../EllipseFit

#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1.csv &> outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2.csv &> outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_38812_39041_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_38812_39041_clean_phase2_Ef.txt


#  Part9_- Example 2. 27036_27279

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1.csv &> outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2.csv &> outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part9_2/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_27036_27279_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part9_2_27036_27279_clean_phase2_Ef.txt






# Part 2 - Example 1. 1847_2021

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
	#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_ViconMrks.txt
	#P2
	#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit

./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_Ef.txt   

./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_Ef.txt   

# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_1847_2021_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_1847_2021_clean_phase2_Ef.txt


./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef.txt   



./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt   

./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt   

./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4598_4797_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4598_4797_clean_phase1_Ef.txt   

./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4598_4797_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4598_4797_clean_phase2_Ef.txt   

# Part 2 - Example 2. 24538_24717


# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2.txt



# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_24538_24717_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_24538_24717_clean_phase2_Ef.txt


# Part 2 - Example 3. 4264_4463


# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2.txt



# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_4264_4463_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_4264_4463_clean_phase2_Ef.txt









#  Part2- Example 4. 843_1042

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_843_1042_clean_phase1.csv &> outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1_843_1042_clean_phase2.csv &> outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase2.txt



# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1042_clean_phase1.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1042_clean_phase2.csv &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1042_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_843_1042_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_843_1042_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_11_24_study_4_Part2_1_843_1042_clean_phase2_Ef.txt










# Part 10 - Example 1. 428_680

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_428_680_clean_phase1.csv &> outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_428_680_clean_phase2.csv &> outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase1.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_428_680_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase2.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_428_680_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_428_680_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_428_680_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_428_680_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_428_680_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_428_680_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_428_680_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_428_680_clean_phase2_Ef.txt




# Part 10 - Example 2. 2132_2486

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1.csv &> outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2.csv &> outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_2132_2486_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_2132_2486_clean_phase2_Ef.txt





#  Part10_1_- Example 3. 6857_7111

# (0). Fit elliptical model to inPts
make ../EllipseFit
#P1
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1.csv &> outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1.txt
#P2
#./../EllipseFit/EllipseFit ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2.csv &> outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2.txt

# (A).Maya 
#		(A.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2.csv &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_ViconMrks.txt

#		(A.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_Ef.txt &> outEfPts/MayaAnimation/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_Ef.txt


# (B).Mathematica 

#		(B.1).inPts vicon
	#P1
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_ViconMrks.txt
	#P2
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2.csv &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_ViconMrks.txt

#		(B.2).inPts Ellipse-Fit
	#P1
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_6857_7111_clean_phase1_Ef.txt
	#P2
#./PathPlanning outEfPts/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_Ef.txt &> outEfPts/MathematicaPlots/2016_12_14_study_4_Part10_1_6857_7111_clean_phase2_Ef.txt




