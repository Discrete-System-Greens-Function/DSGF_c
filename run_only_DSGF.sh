#!/bin/bash

echo "Running DSGF "
./DSGF

if [ $? -eq 139 ]; then
	echo "An error in the simulation occured."
	echo "Verify if the control.txt and the geometry file match the parameters for the simulation."
	echo "If they matched, this simulation may have run out of memory due to number of subvolumes."
	echo "If you are running the code in parallel, disable the option."
	echo "If you are running using the direct option, consider changing to iterative."

fi


