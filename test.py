# run FSM with fixed frequency threshold -f 10000 -d 2 -l 10,000,000
import os

networks = ["string", "YeastPPI", "Rollad"]
network_paths = ["/home/tw83/twang/FSM/networks/string_network_combined_score950.txt","/home/tw83/twang/FSM/networks/YeastPPI", "/home/tw83/twang/FSM/RolladVidal.reformat.uniq.txt"]

for i in range(0,3):
	network = networks[i]
	network_path = network_paths[i]
	for size in range(10, 61, 10):
		our_time_request = "0:15:00"
		our_mem_request = "4G"
		if size <= 30:
			QX_time_request = "0:30:00"
			QX_mem_request = "8G"
		else :
			QX_time_request = "0:30:00"
			QX_mem_request = "32G"
		our_command = "sbatch -p short -t " + our_time_request + " --mem-per-cpu=" + our_mem_request 
		our_command += " --output OurESU." + network + "." + str(size) + ".out " + "--wrap='./OurESU -i " + network_path + " -s " + str(size) + " -d 2 -l 10000000 -f 10000 -u > OurESU." +network+"."+str(size)+".txt' "

		QX_command = "sbatch -p short -t " + QX_time_request + " --mem-per-cpu=" + QX_mem_request 
		QX_command += " --output QX." + network + "." + str(size) + ".out " + "--wrap='./QX -i " + network_path + " -s " + str(size) + "  -l 10000000  -u > QX." +network+"."+str(size)+".txt' "
		os.system(QX_command)
		os.system(our_command)
		break


