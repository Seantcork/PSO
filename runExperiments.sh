
for swarmSize in 16 30 49
	do 
		for topology in gl ri vn ra
			do
				for testFunction in rok ack ras
					do
						for((i = 0; i < 20;i++));
						do
							./pso $topology $swarmSize 10000 $testFunction 30 >> experimentResults.txt 2>&1 &
							wait $!
						done
					done



	done
done
