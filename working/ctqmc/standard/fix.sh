for i in *
do
    echo "current directory:"
    pwd
    cd $i
    echo "job directory:"
    pwd
    sed -i 's///g' atom.config.in
    echo ''
    cd ..
done
