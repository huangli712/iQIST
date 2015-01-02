for i in *
do
    echo "current directory:"
    pwd
    cd $i
    echo "job directory:"
    pwd
    sed -i '' 's/manjushaka/MANJUSHAKA/g' solver.ctqmc.in
    echo ''
    cd ..
done
