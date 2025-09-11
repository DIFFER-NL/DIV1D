
#!/bin/bash


#DESTDIR=/home/unix/derks/codes/div1dbuilds/v304/obj/
PWD=$(pwd)
DIVDIR=$PWD""
CODEDIR=$PWD"/obj"
TESTDIR=$PWD"/int-tests"


cd $TESTDIR
TESTS= */
echo $TESTS

# cd $DIVDIR

echo "run all tests (0), one of them (1-XX), or request input (-1): "
read RUNNER
echo $RUNNER

count=0
for TEST in $(ls -d */);
do
	echo $TEST
	(( count++ ))
	echo $count", run "$TEST":"
	
	if [[ $RUNNER == 0 ]] # run
	then
	echo 'all'
	VAR=1
	
	else
	VAR=0
 		if [[ $RUNNER == -1 ]] # ask
 		then
		echo 'ask'
 		read VAR 
 		else
   			if [[ $RUNNER == $count ]] # single test
   			then
			echo 'single'
			VAR=1
   			fi   
		fi
	fi
	echo "run (1/0): "$VAR
	
	if [[ $VAR == 1 ]]
	then # go into directory to have output.txt and nohup written there
	cd $TEST
	sleep 0.2
	rm div1d_output.txt
	rm -r nohup*.out
        echo $CODEDIR"/div1d.exe < input.txt"
	nohup $CODEDIR"/div1d.exe" < "input.txt" &
	#nohup sh -c 'COUNTER=1; while true; do SIZE=$(stat --printf="%s" nohup.out); if [ "$SIZE" -gt 100000 ]; then cp nohup.out nohup1.out; echo "" > nohup.out; COUNTER=$((COUNTER + 1)); fi; sleep 60; done' >/dev/null 2>&1 &
	cd ..
	fi 
	sleep 0.2
done


sleep 0.2
echo "use the <top> command to see if the div1d.exe processes finish" 


