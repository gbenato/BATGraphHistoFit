#for((ds=3601; ds<=3615; ds++))
#do
#    /home/gbenato/Cuore/Analysis/BATGraphHistoFit/bin/BAT_Reso \
#	--input-path /home/gbenato/Cuore/Data/FakeNature/lineshape_scaling_input/fitresults \
#	--output-path  /home/gbenato/Cuore/Data/FakeNature/lineshape_scaling_output \
#	--label linearreso \
#	--dataset $ds \
#	--precision 3
#done


for((ds=3817; ds<=3818; ds++))
do
    /home/gbenato/Cuore/Analysis/BATGraphHistoFit/bin/BAT_Reso \
	--input-path /home/gbenato/Cuore/Data/TAUP23/lineshape/background/lineshape_fit_results/bkg/plots/fitresults/fitresults \
	--output-path  /home/gbenato/Cuore/Data/TAUP23/lineshape_scaling_output \
	--label linearreso \
	--dataset $ds \
	--precision 3
done
