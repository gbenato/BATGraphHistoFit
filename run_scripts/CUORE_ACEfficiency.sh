#array=(3817 3818 3819 3821 3822 3823 3824 3825 3826 3827)
#array=(3816 3828)
#array=(3817)

#for ds in "${array[@]}"
for((ds=3816; ds<=3828; ds++))
do
    /home/gbenato/Cuore/Analysis/BATGraphHistoFit/bin/BAT_Efficiency \
	--input-path /home/gbenato/Cuore/Data/TAUP23/EfficiencyInput/Eff_ \
	--output-path /home/gbenato/Cuore/Data/TAUP23/ACEfficiency \
	--cfg-path /home/gbenato/Cuore/Analysis/BATGraphHistoFit/cfg/MultiPeaks.cfg \
	--dataset $ds \
	--precision 3 \
	--mode M 
done
