#! /bin/ksh

# Script using GRIB API to extract SML1
# ALBERGEL 2014
####################################################### LOOP on DATES
rm -f H1
rm -f H2
rm -f notGood.txt

for i in `ls /cnrm/vegeo/muciaa/NO_SAVE/US00/ekf_vod10/ISBA_ANALYSIS*_0.nc`
do

echo $i

H1=`ncks -v HO1_1_1 $i | head -11 | tail -1`
H2=`ncks -v HO2_1_1 $i | head -11 | tail -1`

### For LAI+SSM and VOD+SSM
#if [[ $H1 == *"observation WG2 and control variable LAI"* ]] && [[ $H2 == *"observation LAI and control variable LAI"* ]]
### For LAI and VOD
if [[ $H1 == *"observation LAI and control variable LAI"* ]] && [[ $H2 != *"observation WG2 and control variable LAI"* ]]
### For SSM
#if [[ $H1 == *"observation WG2 and control variable LAI"* ]] && [[ $H2 != *"observation WG2 and control variable LAI"* ]]
then
    echo "Good"
else
    echo $i >> notGood.txt
fi
done
