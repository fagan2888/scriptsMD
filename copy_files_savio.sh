#!/bin/bash

### script to automate making folders for nbands/encut convergence tests

while getopts 'neim' flag;
do

        if [ "${flag}" == "m" ];
        then
                T=1010479
                E=6000
                N=2800
                sigma=$(echo "scale=20;0.00008617330337217213" | bc -l)
                sigma=$(echo "scale=16;$sigma * $T" | bc -l)
                cp -p ~/templates/KPOINTS .
                #CHANGE
                cp -p ~/poscars/POSCAR_2fold_15atom POSCAR
                cp -p ~/potcars/POTCAR_MgSiO_LDA_12e_new POTCAR
                sed -e 's/job-name=/job-name=MD'$T'/g' ~/templates/run-savio.job_MD > run-savio.job
                sed -e 's/TEBEG/TEBEG = '$T'/g' -e 's/SIGMA/SIGMA = '$sigma'/g' \
                        -e 's/ENCUT/ENCUT = '$E'/g' \
                        -e 's/NBANDS/NBANDS = '$N'/g' ~/templates/INCAR_MD > INCAR
                sbatch run-savio.job

        elif [ "${flag}" == "n" ];
        then
                name=$(basename $PWD)
                T=2020958
                #T=1347305
                #T=1010479
                E=5000
                for N in 600;
                do
                        sigma=$(echo "scale=20;0.00008617330337217213" | bc -l)
                        sigma=$(echo "scale=16;$sigma * $T" | bc -l)
                        if [ $N -le 99 ];
                        then
                                mkdir NBANDS00$N
                                cd NBANDS00$N
                        elif [ $N -le 999 ];
                        then
                                mkdir NBANDS0$N
                                cd NBANDS0$N
                        else
                                mkdir NBANDS$N
                                cd NBANDS$N
                        fi
                        cp -p ~/templates/KPOINTS .
                        #CHANGE
                        cp -p ~/poscars/POSCAR_12fold POSCAR
                        cp -p ~/potcars/POTCAR_MgSiO_GGA POTCAR
                        sed -e 's/job-name=/job-name='$name'N'$N'/g' ~/templates/run-savio.job > run-savio.job
                        sed -e 's/TEBEG/TEBEG = '$T'/g' \
                                -e 's/SIGMA/SIGMA = '$sigma'/g' \
                                -e 's/ENCUT/ENCUT = '$E'/g' \
                                -e 's/NBANDS/NBANDS = '$N'/g' ~/templates/INCAR > INCAR
                        sbatch run-savio.job
                        cd ..
                done

        elif [ "${flag}" == "e" ];
        then
                T=100000
                N=100
                for E in {3000..5000..1000};
                do
                        sigma=$(echo "scale=20;0.00008617330337217213" | bc -l)
                        sigma=$(echo "scale=16;$sigma * $T" | bc -l)
                        if [ $E -le 99 ];
                        then
                                mkdir ENCUT00$E
                                cd ENCUT00$E
                        elif [ $E -le 999 ];
                        then
                                mkdir ENCUT0$E
                                cd ENCUT0$E
                        else
                                mkdir ENCUT$E
                                cd ENCUT$E
                        fi
                        cp -p ~/templates/KPOINTS .
                        #CHANGE
                        cp -p ~/poscars/POSCAR_2fold POSCAR
                        cp -p ~/potcars/POTCAR_MgSiO_LDA_12e_new POTCAR
                        sed -e 's/job-name=/job-name=E'$E'/g' ~/templates/run-savio.job > run-savio.job
                        sed -e 's/TEBEG/TEBEG = '$T'/g' -e 's/SIGMA/SIGMA = '$sigma'/g' \
                                -e 's/ENCUT/ENCUT = '$E'/g' \
                                -e 's/NBANDS/NBANDS = '$N'/g' ~/templates/INCAR > INCAR
        #                sbatch run-savio.job
                        cd ..
                done
        elif [ "${flag}" == "i" ];
        then
                E=4000
                i=7
                NPAR=$(awk '/NPAR/ {print $3}' ~/templates/INCAR)
                N1=200
                N2=2800
                T1=100000
                T2=1010479
                for T in {250000,500000,750000,900000};
                do
                        if [ $i -le 9 ];
                        then
                                i=00$i
                        elif [ $i -le 99 ];
                        then
                                i=0$i
                        fi
                        N=$(echo "($N2-$N1)/($T2-$T1) * ($T-$T1) + $N1" | bc -l)
                        N=$(echo "($N + 0.5)/1" | bc)
                        N=$(echo "$N - ($N % $NPAR)" | bc)
                        sigma=$(echo "scale=20;0.00008617330337217213" | bc -l)
                        sigma=$(echo "scale=16;$sigma * $T" | bc -l)
                        if [ ! -d MgSiO3_$i ];
                        then
                                mkdir MgSiO3_$i
                        fi
                        cd MgSiO3_$i
                        i=$(echo "$i + 1" | bc)
                        if [ ! -d GGA ];
                        then
                                mkdir GGA
                        fi
                        cd GGA
                        #CHANGE THIS
                        cp -p ~/poscars/POSCAR_2fold POSCAR
                        cp -p ~/potcars/POTCAR_MgSiO_GGA POTCAR
                        cp -p ~/templates/KPOINTS .
                        sed -e 's/job-name=/job-name=T'$T'/g' ~/templates/run-savio.job > run-savio.job
                        sed -e 's/TEBEG/TEBEG = '$T'/g' -e 's/SIGMA/SIGMA = '$sigma'/g' \
                                -e 's/ENCUT/ENCUT = '$E'/g' \
                                -e 's/NBANDS/NBANDS = '$N'/g' ~/templates/INCAR > INCAR
                        sbatch run-savio.job
                        cd ../..
                done
        else
                echo
                echo "Please pass a flag."
                echo "    -n for NBANDS"
                echo "    -e for ENCUT"
                echo "    -i for MgSiO3_..."
                echo
                exit 1
        fi
done
