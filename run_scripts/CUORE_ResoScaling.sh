#array=(3816 3828)

#for ds in "${array[@]}"
#for((ds=3816; ds<=3828; ds++))
for((ds=3601; ds<=3615; ds++))
do
    /home/gbenato/Cuore/Analysis/BATGraphHistoFit/bin/BAT_Reso \
	--input-path /home/gbenato/Cuore/Data/FakeNature/lineshape_scaling_input/fitresults \
	--output-path  /home/gbenato/Cuore/Data/FakeNature/lineshape_scaling_output \
	--label linearreso \
	--dataset $ds \
	--precision 3
done
