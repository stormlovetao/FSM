for size in $(seq 10 10 80);do
	./OurESU -i ./networks/string_network_combined_score950.txt -s $size -d 4 -l 10000000 -f 5000 -u > OurESU.string.${size}.txt
done
