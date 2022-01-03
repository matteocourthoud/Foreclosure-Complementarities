# Make script executable
# chmod a+x Dropbox/Projects/Foreclosure/copyfigures.sh

printf "\nMoving figures...\n\n"

# Move to directory
cd Dropbox/Projects/Foreclosure-complementarities/

# Monopoly timelines
#cp output/alluvial/game_a30g0s70_baseline.png figures/alluvial_monopoly.png
#cp output/timelines/baseline_a30g0s70_welfare.png figures/timeline_monopoly_welfare.png
#cp output/timelines/baseline_a30g0s70_surplus.png figures/timeline_monopoly_surplus.png
#cp output/timelines/baseline_a30g0s70_profits.png figures/timeline_monopoly_profits.png

# Competitive timelines
#cp output/alluvial/game_a70g0s30_baseline.png figures/alluvial_competitive.png
#cp output/timelines/baseline_a70g0s30_welfare.png figures/timeline_competitive_welfare.png
#cp output/timelines/baseline_a70g0s30_surplus.png figures/timeline_competitive_surplus.png
#cp output/timelines/baseline_a70g0s30_profits.png figures/timeline_competitive_profits.png

# Baseline compstats
for value in "margin" "merger" "exit" "profitsLR" "surplusLR" "mpoly"; do
    cp output/compstats/baseline/${value}.png figures/compstats_${value}.png
done

# Effect of learning and bundling
for value in "margin" "merger" "exit" "welfare" "mpoly" "mpolyA" "mpolyB"; do
    for p in "mergers" "learning" "bundling"; do
        #cp output/compstats/no${p}/${value}.png figures/no${p}_${value}.png
        cp output/compstats/no${p}/diff_${value}.png figures/no${p}_${value}_diff.png
    done
done

# Predatory incentives
for value in "welfare" "mpoly" "mpolyA" "mpolyB"; do
    for pred in "bundling" "pricing"; do
        for e in "entry" "exit"; do
            #cp output/compstats/nopred${e}${pred}/${value}.png figures/nopred${pred}${e}_${value}.png
            cp output/compstats/nopred${e}${pred}/diff_${value}.png figures/nopred${pred}${e}_${value}_diff.png
        done
    done
done

# Policy
for value in "margin" "merger" "exit" "welfare" "mpoly"; do
    for policy in "limitedmergers" "datasharing" "limitedbundling"; do
        #cp output/compstats/${policy}/${value}.png figures/${policy}_${value}.png
        cp output/compstats/${policy}/diff_${value}.png figures/${policy}_${value}_diff.png
    done
done

# Terminate
exit
