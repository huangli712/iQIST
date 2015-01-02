for i in *
do
    echo "current directory:"
    pwd
    cd $i
    echo "job directory:"
    pwd
    sed -i '' 's/pansy/PANSY/g' atom.config.in
    echo ''
    cd ..
done
