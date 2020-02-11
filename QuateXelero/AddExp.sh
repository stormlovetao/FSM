#!/bin/bash


mkdir ecoli_5
sbatch -p short -t 2:00:00 -e AddExp.run5.elog -o AddExp.run5.olog --mem-per-cpu=20G --wrap='./QuateXelero -i networks/ecoli -s 5 -r 1000 -o ./ecoli_5/ > ./ecoli_5/log.txt'

mkdir ecoli_6
sbatch -p short -t 2:00:00 -e AddExp.run6.elog -o AddExp.run6.olog --mem-per-cpu=20G --wrap='./QuateXelero -i networks/ecoli -s 6 -r 1000 -o ./ecoli_6/ > ./ecoli_6/log.txt'


mkdir ecoli_7
sbatch -p short -t 2:00:00 -e AddExp.run7.elog -o AddExp.run7.olog --mem-per-cpu=30G --wrap='./QuateXelero -i networks/ecoli -s 7 -r 1000 -o ./ecoli_7/ > ./ecoli_7/log.txt'

mkdir ecoli_8
sbatch -p short -t 4:00:00 -e AddExp.run8.elog -o AddExp.run8.olog --mem-per-cpu=40G --wrap='./QuateXelero -i networks/ecoli -s 8 -r 1000 -o ./ecoli_8/ > ./ecoli_8/log.txt'


