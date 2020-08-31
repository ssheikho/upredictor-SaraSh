#!/bin/bash  
# int-or-string.sh

make

declare -i reachCount
# The script will treat subsequent occurrences of "number" as an integer.	
# ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_515_717_clean_phase1.csv&>out.txt
# ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part2_1/2016_11_24_study_4_Part2_1_40174_40395_clean_phase1.csv&>out2.txt
# ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_11_24_study_4_Part3_1/2016_11_24_study_4_Part3_1_253_479_clean_phase1.csv&>out3.txt
# ./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_36294_36593_clean_phase1.csv&>out10.txt
 
#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part12_1/2016_12_14_study_4_Part12_1_34559_34796_clean_phase1.csv&>out12_7.txt

./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase1.csv &>out10-1.txt


#./PathPlanning ~/Dropbox/Reaching\ Study/Gilwoo_Automated_v2/csv_only_manual_inspection/2016_12_14_study_4_Part10_1/2016_12_14_study_4_Part10_1_428_680_clean_phase2.csv &>out10-2.txt

