# bahllab
```
interact --ntasks 8 --mem 64G
NOTEBOOKPORT=8888
IPUSED=$(hostname -i)
echo "NOTEBOOKPORT is " $NOTEBOOKPORT
echo "IPUSED is " $IPUSED
module load Anaconda3/2020.02
jupyter-lab --port $NOTEBOOKPORT --ip=$IPUSED --no-browser
```

```
ssh -N -L 8888:10.2.1.91:8888 gs69042@sapelo2.gacrc.uga.edu
```
