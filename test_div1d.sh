#DESTDIR=/home/unix/derks/codes/div1dbuilds/v304/obj/
PWD=$(pwd)
DIVDIR=$PWD""
CODEDIR=$PWD"/obj"
TESTDIR=$PWD"/unit-tests"

# echo "PWD before changing directory"
#echo $PWD
#cd $CODEDIR  # directory change
# echo "PWD after changing directory"
#echo $PWD
# other directories only depend on PWD at their definition
#echo $DIVDIR
#echo $CODEDIR
#echo $TESTDIR

echo "The following tests are only to test if runs start and finish "
echo "Finished tests CANNOT be used to judge the correctness of outcomes "

echo "1 test regular input files DIV1D: (1/0)"
read VAR 
if [[ $VAR == 1 ]]
then
# go into directory to have output.txt and nohup written there
cd $TESTDIR"/test1"
nohup $CODEDIR"/div1d.exe" < $TESTDIR"/test1/input.txt" &
cd $DIVDIR
fi 
sleep 0.2
 
echo  "2 test dynamic input files DIV1D: (1/0)"
read VAR 
if [[ $VAR == 1 ]]
then
cd $TESTDIR"/test2"
nohup $CODEDIR"/div1d.exe" < $TESTDIR"/test2/input.txt" &
cd $DIVDIR
fi
sleep 0.2


echo  "3 test external call from C: (1/0)"
read VAR 
if [[ $VAR == 1 ]]
then
cd $TESTDIR"/testlib"
echo  "testing external library call from C "
nohup $CODEDIR"/test_libdiv1d.exe" &
cd $DIVDIR
fi
sleep 0.2
echo "use the <top> command to see if the div1d.exe processes finish" 
