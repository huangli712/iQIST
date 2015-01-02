for i in *
do
    echo "current directory:"
    pwd
    cd $i
    echo "job directory:"
    pwd
    sed -i '' 's/begonia, lavender, pansy and manjushaka/BEGONIA, LAVENDER, PANSY, and MANJUSHAKA/g' atom.config.in
    echo ''
    cd ..
done
