#!/bin/bash
for size in $(seq 10 10 80);do
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo OurESU.string.${size}.out "./OurESU -i ./networks/string_network_combined_score950.txt -s $size -d 4 -l 10000000 -f 5000 -u > OurESU.string.${size}.txt"
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo QX.string.${size}.out "./QX -i ./networks/string_network_combined_score950.txt -s $size -l 10000000 -u > QX.string.${size}.txt"
done

for size in $(seq 10 10 80);do
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo OurESU.YeastPPI.${size}.out "./OurESU -i ./networks/YeastPPI -s $size -d 2 -l 10000000 -f 5000 -u > OurESU.YeastPPI.${size}.txt"
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo QX.YeastPPI.${size}.out "./QX -i ./networks/YeastPPI -s $size -l 10000000 -u > QX.YeastPPI.${size}.txt"
done

for size in $(seq 10 10 80);do
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo OurESU.yeast.${size}.out "./OurESU -i ./networks/yeast -s $size -d 2 -l 10000000 -f 5000 > OurESU.yeast.${size}.txt"
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo QX.yeast.${size}.out "./QX -i ./networks/yeast -s $size -l 10000000 > QX.yeast.${size}.txt"
done

for size in $(seq 10 10 80);do
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo OurESU.Rollad.${size}.out "./OurESU -i RolladVidal.reformat.uniq.txt -s $size -d 4 -l 10000000 -f 5000 -u > OurESU.Rollad.${size}.txt"
	bsub -q short -W 4:00 -M 30000 -R 'rusage[mem=30000]' -oo QX.Rollad.${size}.out "./QX -i RolladVidal.reformat.uniq.txt -s $size -l 10000000  -u > QX.Rollad.${size}.txt"
done
