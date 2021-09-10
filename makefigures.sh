# Make script executable
# chmod a+x Dropbox/Projects/Foreclosure/copyfigures.sh

printf "\nMoving figures...\n\n"

# Move to directory
cd Dropbox/Projects/Foreclosure-complementarities/

# Baseline
for value in "M" "PI" "CS" "W" "P" "Mg" "B" "E"; do
    cp output/figures/compstats/baseline/${value}_lbd.png figures/compstats_${value}.png
done

# Monopoly
cp output/figures/alluvial/lbd_sigma.15_rho.25.png figures/alluvial_monopoly.png
cp output/data/statestats/lbd_sigma.15_rho.25.csv figures/statestats_monopoly.csv
cp output/figures/timelines/lbd_sigma.15_rho.25_W.png figures/timeline_monopoly_W.png
cp output/figures/timelines/lbd_sigma.15_rho.25_CS.png figures/timeline_monopoly_CS.png
cp output/figures/timelines/lbd_sigma.15_rho.25_PI.png figures/timeline_monopoly_PI.png

# Competitive
cp output/figures/alluvial/lbd_sigma.25_rho.75.png figures/alluvial_competitive.png
cp output/data/statestats/lbd_sigma.25_rho.75.csv figures/statestats_competitive.csv
cp output/figures/timelines/lbd_sigma.25_rho.75_W.png figures/timeline_competitive_W.png
cp output/figures/timelines/lbd_sigma.25_rho.75_CS.png figures/timeline_competitive_CS.png
cp output/figures/timelines/lbd_sigma.25_rho.75_PI.png figures/timeline_competitive_PI.png

# Effect of learning and bundling
for value in "P" "B" "E" "M" "CS" "W"; do
    for p in "bundling" "learning"; do
        cp output/figures/compstats/${p}/${value}_lbd.png figures/no${p}_${value}_diff.png
        cp output/figures/compstats/baseline/${value}_lbd_${p}.n.png figures/no${p}_${value}.png
    done
done


# Predatory incentives
for value in "P" "B" "E" "M" "CS" "W"; do
    for pred in "predbundl" "predprice"; do
        for e in "entry" "exit"; do
            cp output/figures/compstats/${pred}\ ${e}/${value}_lbd.png figures/no${pred}_${e}_${value}_diff.png
            cp output/figures/compstats/baseline/${value}_lbd_${e}.no${pred}.png figures/no${pred}_${e}_${value}.png
        done
    done
done

# Policy
for value in "P" "B" "E" "M" "CS" "W"; do
    for policy in "limitmergers" "limitbundling" "datasharing"; do
        cp output/figures/compstats/${policy}/${value}_lbd.png figures/${policy}_${value}_diff.png
    done
    cp output/figures/compstats/baseline/${value}_lbd_mergers.limited.png figures/limitmergers_${value}.png
    cp output/figures/compstats/baseline/${value}_lbd_bundling.limited.png figures/limitbundling_${value}.png
    cp output/figures/compstats/baseline/${value}_lbd_maxexpgap.100.png figures/datasharing_${value}.png
done

# Terminate
exit
