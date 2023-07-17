#!/bin/bash

cmd1='CHKERRQ(IERR)'
cmd2='if (IERR .ne. 0) then;call PetscErrorF(IERR);return;endif'
files='mod_petsc.F mod_semi_implicit.F mod_non_hydro.F mod_wavesetup.F mod_action_im.F'

echo "========================================"
echo " Convert the PetSc CHKERRQ:"
echo "    --- 1  $cmd1 -> $cmd2"
echo "    --- 2  $cmd2 -> $cmd1"
read choice

if [ "$choice" == "1" ]; then
    for file in $files; do
      echo $file
      sed -i "s/${cmd1}/${cmd2}/g" $file
    done
elif [ "$choice" == "2" ]; then
    for file in $files; do
      echo $file
      sed -i "s/${cmd2}/${cmd1}/g" $file
    done
else
    # Handle invalid input
    echo "Error: Invalid choice. Please enter either 1 or 2."
fi


